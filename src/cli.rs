use std::path::PathBuf;
use structopt::StructOpt;
use strum_macros::EnumString;

// This module contains the `Command` enum, which is the wrapper for all subcommands of the `rpt` binary.

#[derive(StructOpt)]
#[structopt(
    about = "A set of command line utilities for phylogeny related bioinformatics tasks.",
    author = "Ragnar Groot Koerkamp",
    name = "Rust-Phylogeny-Tools"
)]
pub(crate) struct Rpt {
    //#[structopt(long, short, help = "Verbose output.")]
    //pub(crate) verbose: bool,
    #[structopt(subcommand)]
    pub(crate) cmd: Command,
}

#[derive(EnumString)]
pub(crate) enum PhylogenyMethod {
    UPGMA,
    NeighborJoining,
}

#[derive(StructOpt)]
pub(crate) enum Command {
    // Reconstruct the phylogeny, given a distance matrix
    Phylogeny {
        #[structopt(long = "method")]
        method: PhylogenyMethod,

        // Input .phylip distance matrix.
        #[structopt(parse(from_os_str))]
        input: PathBuf,

        // Path to store the phylogeny in .newick format, or stdout otherwise.
        #[structopt(parse(from_os_str))]
        output: Option<PathBuf>,
    },

    // Given two .newick files, compute the RF-distance between them and print it to stdout.
    RobinsonFoulds {
        #[structopt(parse(from_os_str))]
        newick_1: PathBuf,
        #[structopt(parse(from_os_str))]
        newick_2: PathBuf,

        // Path to store the distance, or stdout otherwise.
        #[structopt(parse(from_os_str))]
        output: Option<PathBuf>,
    },
}
