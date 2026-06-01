use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// A hierarchical linkage matrix: `n-1` agglomeration rows in scipy order, each
/// `[left, right, height, size]`. Cluster ids `< n` are original observations,
/// ids `>= n` refer to the cluster formed at row `id - n`.
pub struct Linkage {
    pub left: Vec<usize>,
    pub right: Vec<usize>,
    pub height: Vec<f64>,
    pub n: usize,
}

impl Linkage {
    #[must_use]
    pub fn rows(&self) -> usize {
        self.n - 1
    }
}

/// Read the linkage matrix written by rsomics-upgma (and scipy): tab-separated
/// `left right height size`, `n-1` rows, no header.
pub fn read_linkage(input: &Path) -> Result<Linkage> {
    let reader: Box<dyn BufRead> = if input.as_os_str() == "-" {
        Box::new(BufReader::new(std::io::stdin().lock()))
    } else {
        let f = std::fs::File::open(input)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
        Box::new(BufReader::new(f))
    };

    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut height = Vec::new();

    for (r, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() != 4 {
            return Err(RsomicsError::InvalidInput(format!(
                "linkage row {r} has {} columns, expected 4 (left right height size)",
                cells.len()
            )));
        }
        let parse_idx = |s: &str| -> Result<usize> {
            s.trim()
                .parse::<f64>()
                .map(|v| v as usize)
                .map_err(|e| RsomicsError::InvalidInput(format!("row {r} bad cluster id: {e}")))
        };
        left.push(parse_idx(cells[0])?);
        right.push(parse_idx(cells[1])?);
        height.push(
            cells[2]
                .trim()
                .parse::<f64>()
                .map_err(|e| RsomicsError::InvalidInput(format!("row {r} bad height: {e}")))?,
        );
    }

    let rows = height.len();
    if rows == 0 {
        return Err(RsomicsError::InvalidInput("empty linkage matrix".into()));
    }
    let n = rows + 1;

    Ok(Linkage {
        left,
        right,
        height,
        n,
    })
}
