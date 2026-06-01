use std::io::Write;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub mod cluster;
pub mod linkage;

pub use cluster::{Criterion, fcluster};
pub use linkage::{Linkage, read_linkage};

/// Read a length-`n-1` monocrit statistic vector (one f64 per line).
pub fn read_monocrit(path: &Path, rows: usize) -> Result<Vec<f64>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let v: Vec<f64> = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.trim()
                .parse::<f64>()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad monocrit value: {e}")))
        })
        .collect::<Result<_>>()?;
    if v.len() != rows {
        return Err(RsomicsError::InvalidInput(format!(
            "monocrit has {} values, expected {rows} (n-1)",
            v.len()
        )));
    }
    Ok(v)
}

/// Write one cluster label per line, observation order 0..n. Matches scipy's
/// `T` array ordering.
pub fn write_labels(labels: &[i32], out: &mut dyn Write) -> Result<()> {
    for c in labels {
        writeln!(out, "{c}").map_err(RsomicsError::Io)?;
    }
    Ok(())
}
