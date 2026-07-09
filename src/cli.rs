use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, HelpSpec, Origin};

use rsomics_fcluster::{Criterion, fcluster, read_linkage, read_monocrit, write_labels};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum CriterionArg {
    Distance,
    Maxclust,
    Inconsistent,
    Monocrit,
}

impl From<CriterionArg> for Criterion {
    fn from(c: CriterionArg) -> Self {
        match c {
            CriterionArg::Distance => Criterion::Distance,
            CriterionArg::Maxclust => Criterion::MaxClust,
            CriterionArg::Inconsistent => Criterion::Inconsistent,
            CriterionArg::Monocrit => Criterion::Monocrit,
        }
    }
}

#[derive(Parser, Debug)]
#[command(name = "rsomics-fcluster", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Linkage matrix TSV (left right height size), n-1 rows. "-" for stdin.
    pub input: PathBuf,

    #[arg(long = "threshold")]
    threshold: f64,

    #[arg(long = "criterion", value_enum, default_value = "distance")]
    criterion: CriterionArg,

    #[arg(long = "depth", default_value = "2")]
    depth: usize,

    #[arg(long = "monocrit")]
    monocrit: Option<PathBuf>,

    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let z = read_linkage(&self.input)?;
        let crit: Criterion = self.criterion.into();

        let monocrit = if crit == Criterion::Monocrit {
            let path = self.monocrit.as_ref().ok_or_else(|| {
                RsomicsError::InvalidInput(
                    "--criterion monocrit requires --monocrit FILE (n-1 statistics)".into(),
                )
            })?;
            Some(read_monocrit(path, z.rows())?)
        } else {
            None
        };

        let labels = fcluster(&z, self.threshold, crit, self.depth, monocrit.as_deref());

        let mut out: Box<dyn std::io::Write> = if self.output == "-" && self.common.json {
            Box::new(std::io::sink())
        } else if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        write_labels(&labels, &mut out)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Form flat clusters from a hierarchical linkage matrix by cutting it.",
    origin: Some(Origin {
        upstream: "scipy.cluster.hierarchy.fcluster",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &[
        "<linkage.tsv> --threshold T [--criterion distance|maxclust|inconsistent|monocrit] [-o labels.tsv]",
    ],
    sections: &[],
    examples: &[
        Example {
            description: "Cut at cophenetic distance 1.1",
            command: "rsomics-fcluster Z.tsv --threshold 1.1 --criterion distance",
        },
        Example {
            description: "Form at most 4 clusters",
            command: "rsomics-fcluster Z.tsv --threshold 4 --criterion maxclust",
        },
        Example {
            description: "Cut by inconsistency coefficient with depth 3",
            command: "rsomics-fcluster Z.tsv --threshold 0.9 --criterion inconsistent --depth 3",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
