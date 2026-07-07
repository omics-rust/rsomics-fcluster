use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fcluster"))
}

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run(args: &[&str]) -> String {
    let out = Command::new(bin())
        .args(args)
        .output()
        .expect("spawn rsomics-fcluster");
    assert!(
        out.status.success(),
        "rsomics-fcluster {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8 stdout")
}

/// OUR labels must equal scipy's committed labels, byte-for-byte. Always runs.
#[test]
fn golden_labels_match_scipy() {
    let z = golden("small_linkage.tsv");
    let zs = z.to_str().unwrap();
    let cases: &[(&str, &[&str])] = &[
        (
            "small_distance_t1.5.expected",
            &["--threshold", "1.5", "--criterion", "distance"],
        ),
        (
            "small_distance_t3.0.expected",
            &["--threshold", "3.0", "--criterion", "distance"],
        ),
        (
            "small_maxclust_4.expected",
            &["--threshold", "4", "--criterion", "maxclust"],
        ),
        (
            "small_maxclust_6.expected",
            &["--threshold", "6", "--criterion", "maxclust"],
        ),
        // Degenerate maxclust: t=0 requests at most zero clusters. scipy returns
        // every point in its own cluster, labelled in DFS-reach order.
        (
            "small_maxclust_0.expected",
            &["--threshold", "0", "--criterion", "maxclust"],
        ),
        (
            "small_inconsistent_t1.2_d2.expected",
            &[
                "--threshold",
                "1.2",
                "--criterion",
                "inconsistent",
                "--depth",
                "2",
            ],
        ),
        (
            "small_inconsistent_t0.8_d3.expected",
            &[
                "--threshold",
                "0.8",
                "--criterion",
                "inconsistent",
                "--depth",
                "3",
            ],
        ),
    ];
    for (expected, flags) in cases {
        let mut args = vec![zs];
        args.extend_from_slice(flags);
        let got = run(&args);
        let want = std::fs::read_to_string(golden(expected)).unwrap();
        assert_eq!(got, want, "mismatch for {expected} with {flags:?}");
    }
}

/// Live differential against scipy if a python with scipy is reachable. Loud-skip
/// when absent so CI (no scipy) stays green on the committed golden above.
#[test]
fn live_scipy_differential() {
    let Some(py) = find_python_with_scipy() else {
        eprintln!("SKIP live_scipy_differential: no python with scipy on PATH / SCIPY_PYTHON");
        return;
    };

    let scratch = std::env::temp_dir().join("rsomics-fcluster-compat");
    std::fs::create_dir_all(&scratch).unwrap();
    let z_path = scratch.join("linkage.tsv");

    let gen_src = format!(
        r#"
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
rng = np.random.default_rng(777)
pts = np.vstack([rng.normal(c, 0.7, (15,4)) for c in ([0,0,0,0],[5,0,0,0],[0,5,0,0],[0,0,5,0],[0,0,0,5])])
Z = linkage(pdist(pts), method="average")
with open(r"{z}", "w") as f:
    for row in Z:
        f.write(f"{{int(row[0])}}\t{{int(row[1])}}\t{{float(row[2]):.17f}}\t{{int(row[3])}}\n")
import json
out = {{}}
for name, kw in [
    ("d2.0", dict(t=2.0, criterion="distance")),
    ("d4.0", dict(t=4.0, criterion="distance")),
    ("mc5", dict(t=5, criterion="maxclust")),
    ("mc8", dict(t=8, criterion="maxclust")),
    ("in1.0d2", dict(t=1.0, criterion="inconsistent", depth=2)),
    ("in0.6d4", dict(t=0.6, criterion="inconsistent", depth=4)),
]:
    out[name] = fcluster(Z, **kw).tolist()
print(json.dumps(out))
"#,
        z = z_path.display()
    );

    let scipy_out = Command::new(&py)
        .arg("-c")
        .arg(&gen_src)
        .output()
        .expect("spawn scipy");
    assert!(
        scipy_out.status.success(),
        "scipy gen failed: {}",
        String::from_utf8_lossy(&scipy_out.stderr)
    );
    let json = String::from_utf8(scipy_out.stdout).unwrap();
    let cases: &[(&str, &[&str])] = &[
        ("d2.0", &["--threshold", "2.0", "--criterion", "distance"]),
        ("d4.0", &["--threshold", "4.0", "--criterion", "distance"]),
        ("mc5", &["--threshold", "5", "--criterion", "maxclust"]),
        ("mc8", &["--threshold", "8", "--criterion", "maxclust"]),
        (
            "in1.0d2",
            &[
                "--threshold",
                "1.0",
                "--criterion",
                "inconsistent",
                "--depth",
                "2",
            ],
        ),
        (
            "in0.6d4",
            &[
                "--threshold",
                "0.6",
                "--criterion",
                "inconsistent",
                "--depth",
                "4",
            ],
        ),
    ];
    let zs = z_path.to_str().unwrap();
    for (key, flags) in cases {
        let mut args = vec![zs];
        args.extend_from_slice(flags);
        let got: Vec<i32> = run(&args).lines().map(|l| l.parse().unwrap()).collect();
        let want = scipy_labels(&json, key);
        assert_eq!(got, want, "live scipy mismatch for {key} {flags:?}");
    }
}

fn scipy_labels(json: &str, key: &str) -> Vec<i32> {
    let needle = format!("\"{key}\": [");
    let start = json.find(&needle).unwrap() + needle.len();
    let end = json[start..].find(']').unwrap() + start;
    json[start..end]
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect()
}

fn find_python_with_scipy() -> Option<String> {
    let mut candidates: Vec<String> = Vec::new();
    if let Ok(p) = std::env::var("SCIPY_PYTHON") {
        candidates.push(p);
    }
    candidates.push("python3".into());
    candidates.push("python".into());
    for c in candidates {
        let ok = Command::new(&c)
            .args(["-c", "import scipy"])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if ok {
            return Some(c);
        }
    }
    None
}
