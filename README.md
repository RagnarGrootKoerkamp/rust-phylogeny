# Rust-Phylogeny, a Rust library for phylogenetic trees

[![Crates.io](https://img.shields.io/crates/d/phylogeny.svg)](https://crates.io/crates/phylogeny)
[![Crates.io](https://img.shields.io/crates/v/phylogeny.svg)](https://crates.io/crates/phylogeny)
[![Crates.io](https://img.shields.io/crates/l/phylogeny.svg)](https://crates.io/crates/phylogeny)
[![GitHub Workflow Status (branch)](https://img.shields.io/github/workflow/status/RagnarGrootKoerkamp/rust-phylogeny/CI/master?label=tests)](https://github.com/ragnargrootkoerkamp/rust-phylogeny/actions)

Currently this library provides three algorithms:

- [UPGMA](https://en.wikipedia.org/wiki/UPGMA) and
  [Nieghbor-Joining](https://en.wikipedia.org/wiki/Neighbor_joining) for
  constructing a phylogeny from a distance matrix.
- [Robinson-Foulds
  metric](https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric) for
  computing a distance between to phylogenetic trees.

Documentation is at https://docs.rs/phylogeny.

This package is inspired by [rust-bio](https://github.com/rust-bio/rust-bio), which you
should also check out.

**This package is not stable.**

## License

Licensed under the MIT license http://opensource.org/licenses/MIT. This project may not be copied, modified, or distributed except according to those terms.
