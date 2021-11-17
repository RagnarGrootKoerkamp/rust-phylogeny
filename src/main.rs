#![warn(missing_docs)]

//! This provides the `rpt` binary, which wraps the algorithms in the `phylogeny` crate.
//!
//! Run `rpt -h` for help.
//!
//! Currently available subcommands are:
//! - `rpt phylogeny --method {upgma,neighbor-joining} <path/to/file.phylip> <path/to/output.newick>`
//! - `rpt robinson-foulds <path/to/input_1.newick> <path/to/input_2.newick>`

use std::fs::File;

use cli::Command;
use phylogeny::*;
use std::io::prelude::*;
use structopt::StructOpt;

mod cli;

fn main() -> phylogeny::Result<()> {
    let args = cli::Rpt::from_args();
    match args.cmd {
        Command::Phylogeny {
            method,
            input,
            output,
        } => {
            let dm = DistanceMatrix::from_file(&input)?;
            let phylogeny = (match method {
                cli::PhylogenyMethod::UPGMA => upgma,
                cli::PhylogenyMethod::NeighborJoining => neighbor_joining,
            })(&dm);
            match output {
                Some(path) => {
                    let mut f = File::create(path)?;
                    f.write_all(phylogeny.to_string().as_bytes())?;
                    f.write(b"\n")?;
                }
                None => {
                    println!("{}", phylogeny.to_string());
                }
            }
        }
        Command::RobinsonFoulds {
            ref newick_1,
            ref newick_2,
            output,
        } => {
            let dist = robinson_foulds_distance(
                &Phylogeny::from_file(newick_1)?,
                &Phylogeny::from_file(newick_2)?,
            );
            match output {
                Some(path) => {
                    let mut f = File::create(path)?;
                    f.write_all(dist.to_string().as_bytes())?;
                    f.write(b"\n")?;
                }
                None => {
                    println!("{}", dist.to_string());
                }
            }
        }
    }
    Ok(())
}
