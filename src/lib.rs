#![warn(missing_docs)]
#![feature(destructuring_assignment)]

//! This library provides algorithms and helper functions for constructing and comparing phylogenies.
//!
//! Currently, it contains these datatypes:
//! - DistanceMatrix, representing a .phylip square symmetric distance matrix.
//! - Phylogeny, representing a .newick file phylogeny.
//!
//! And it contains these algorithms:
//! - UPGMA, a phylogeny reconstruction method.
//! - Neighbor-Joining, another phylogeny reconstruction method.
//! - Robinson-Foulds metrid, a distance metric between phylogenies.
//!
//! There is an accompanying binary `rpt` to call these algorithms from the command line.
//! Run `cargo install phylogeny` followed by `rpt --help` for more details.

use itertools::zip;
use itertools::Itertools;
use ordered_float::NotNan;
use std::cmp::max;
use std::cmp::min;
use std::cmp::Reverse;
use std::collections::BTreeSet;
use std::collections::BinaryHeap;
use std::collections::HashSet;
use std::ffi::OsStr;
use std::fs;
use std::path::Path;
use std::str::FromStr;

type Error = Box<dyn std::error::Error>;

/// A convenience result type.
pub type Result<T> = std::result::Result<T, Error>;

/// A symmetric square matrix representing a .phylip file.
pub struct DistanceMatrix {
    names: Vec<String>,
    distances: Vec<Vec<f32>>,
}

impl DistanceMatrix {
    /// The number of sequences (rows) in the matrix.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Create a new DistanceMatrix.
    ///
    /// `distances` must be a symmetric square matrix with the same number of rows as `names`.
    pub fn new(names: Vec<String>, distances: Vec<Vec<f32>>) -> Self {
        Self { names, distances }
    }

    /// Read a `.phylip` file into a DistanceMatrix.
    pub fn from_file(p: &Path) -> Result<Self> {
        assert!(p.extension() == Some(OsStr::new("phylip")));
        fs::read_to_string(p)?.parse()
    }

    /// Write the DistanceMatrix to a `.phylip` file.
    pub fn to_file(self: &Self, p: &Path) -> Result<()> {
        assert!(p.extension() == Some(OsStr::new("phylip")));
        Ok(fs::write(p, self.to_string())?)
    }
}

impl FromStr for DistanceMatrix {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let mut lines = s.lines();
        let n = lines.next().ok_or("Expected n")?.parse()?;

        let mut d = DistanceMatrix::new(Vec::with_capacity(n), Vec::with_capacity(n));

        for line in lines.into_iter() {
            let mut it = line.split_ascii_whitespace();
            d.names.push(it.next().unwrap().to_string());
            d.distances
                .push(it.map(|chars| chars.parse().unwrap()).collect());
        }

        Ok(d)
    }
}

impl ToString for DistanceMatrix {
    fn to_string(&self) -> String {
        let mut s = String::new();
        s += &self.len().to_string();
        s += "\n";
        for (name, ds) in zip(&self.names, &self.distances) {
            s += &name;
            for d in ds {
                s += " ";
                s += &d.to_string();
            }
            s += "\n";
        }
        s
    }
}

/// A phylogeny representing a .newick file.
#[derive(PartialEq, Debug)]
pub struct Phylogeny {
    /// The name of the current node.
    ///
    /// Can be empty for internal nodes.
    name: String,

    /// The children of the current node.
    ///
    /// Empty for leafs, and distances to the parent are optional.
    children: Vec<(Phylogeny, Option<f32>)>,
}

impl Phylogeny {
    /// Create a new phylogeny node from a name and list of children.
    pub fn new(name: &str, children: Vec<(Phylogeny, Option<f32>)>) -> Phylogeny {
        Phylogeny {
            name: name.to_string(),
            children,
        }
    }

    /// Create a new leaf node.
    pub fn new_leaf(name: &str) -> Phylogeny {
        Phylogeny::new(name, vec![])
    }

    /// Join two trees into a new named parent node, with given edge lengths.
    pub fn join_with_name(name: &str, l: Phylogeny, dl: f32, r: Phylogeny, dr: f32) -> Phylogeny {
        Phylogeny::new(name, vec![(l, Some(dl)), (r, Some(dr))])
    }

    /// Join two trees into a new anonymous parent node, with given edge lengths.
    pub fn join(l: Phylogeny, dl: f32, r: Phylogeny, dr: f32) -> Phylogeny {
        Phylogeny::join_with_name("", l, dl, r, dr)
    }

    /// Read a `.newick` file into a Phylogeny.
    pub fn from_file(p: &Path) -> Result<Self> {
        assert!(p.extension() == Some(OsStr::new("newick")));
        fs::read_to_string(p)?.parse()
    }

    /// Write the phylogeny into a `.newick` file.
    pub fn to_file(self: &Self, p: &Path) -> Result<()> {
        assert!(p.extension() == Some(OsStr::new("newick")));
        Ok(fs::write(p, self.to_string())?)
    }

    fn to_string_impl(&self) -> String {
        if self.children.is_empty() {
            self.name.to_string()
        } else {
            "(".to_string()
                + &self
                    .children
                    .iter()
                    .map(|child| {
                        let (c, d) = child;
                        if let Some(d) = d {
                            (Self::to_string_impl(c) + ":" + &d.to_string()).to_string()
                        } else {
                            Self::to_string_impl(c).to_string()
                        }
                    })
                    .collect::<Vec<String>>()
                    .join(",")
                + ")"
                + &self.name
        }
    }

    fn from_str_impl(mut s: &str) -> Result<(Self, &str)> {
        let mut children = vec![];
        if let Some(_s) = s.strip_prefix('(') {
            s = _s;
            loop {
                let c;
                (c, s) = Self::from_str_impl(s)?;

                let len = if let Some(_s) = s.strip_prefix(':') {
                    s = _s;
                    // The length of the branch ends with , or )
                    let end = s.find(|c| c == ',' || c == ')').ok_or("No , or ) found")?;
                    let len;
                    (len, s) = s.split_at(end);
                    Some(
                        len.parse()
                            .map_err(|e: std::num::ParseFloatError| -> String {
                                "Could not parse length: ".to_string() + &e.to_string()
                            })?,
                    )
                } else {
                    None
                };

                children.push((c, len));
                let done = s.starts_with(')');
                s = s.split_at(1).1;
                if done {
                    break;
                }
            }
        }

        let end = s
            .find(|c| c == ':' || c == ';' || c == ')' || c == ',')
            .ok_or("No : or ; found")?;
        let (name, s) = s.split_at(end);
        Ok((Phylogeny::new(name, children), s))
    }
}

impl ToString for Phylogeny {
    fn to_string(&self) -> String {
        Phylogeny::to_string_impl(&self) + ";"
    }
}

impl FromStr for Phylogeny {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        Ok(Phylogeny::from_str_impl(s)?.0)
    }
}

/// Run the UPGMA algorithm.
///
/// Compute the rooted phylogeny given a `DistanceMatrix`.
pub fn upgma(distances: &DistanceMatrix) -> Phylogeny {
    let n = distances.len();
    // Tuples of (phylogeny, number of leafs in tree, distance to leaf)
    let mut parts: Vec<Option<(Phylogeny, usize, f32)>> = distances
        .names
        .iter()
        .map(|name| Some((Phylogeny::new_leaf(name), 1, 0.0)))
        .collect();
    let mut distances = distances.distances.clone();
    let mut heap = BinaryHeap::with_capacity(n * n);
    for (i, ds) in distances.iter().enumerate() {
        for (j, &d) in ds.iter().enumerate() {
            if i == j {
                continue;
            }
            heap.push(Reverse((NotNan::new(d).unwrap(), i, j)));
        }
    }

    while let Some(Reverse((d, i, j))) = heap.pop() {
        assert!(i != j);
        if parts[i].is_none() || parts[j].is_none() || distances[i][j] != *d {
            continue;
        }
        let (i, j) = (min(i, j), max(i, j));
        if let (Some((pi, si, di)), Some((pj, sj, dj))) = (parts[i].take(), parts[j].take()) {
            // Merge phylogeny i and j, adding theirs sizes. Store the result in the lower of the two indices.
            parts[i] = Some((
                Phylogeny::join(pi, *d / 2. - di, pj, *d / 2. - dj),
                si + sj,
                *d / 2.,
            ));
            // Update the distances to other nodes
            for k in 0..n {
                if k == i || k == j {
                    continue;
                }
                let dk =
                    (si as f32 * distances[i][k] + sj as f32 * distances[j][k]) / (si + sj) as f32;
                distances[i][k] = dk;
                distances[k][i] = dk;
                heap.push(Reverse((NotNan::new(dk).unwrap(), k, i)));
            }
        }
    }

    parts[0].take().unwrap().0
}

/// Run the Neighbor-Joining algorithm
///
/// Computes an unrooted phylogony given a `DistanceMatrix`.
// The root vertex is arbitrarily chosen as the middle in between the last two remaining vertices.
pub fn neighbor_joining(distances: &DistanceMatrix) -> Phylogeny {
    // Tuples of (phylogeny, number of leafs in tree, distance to leaf delta)
    let mut parts: Vec<Option<Phylogeny>> = distances
        .names
        .iter()
        .map(|name| Some(Phylogeny::new_leaf(name)))
        .collect();
    let mut distances = distances.distances.clone();

    let mut active: Vec<usize> = (0..distances.len()).collect();

    while active.len() > 2 {
        let sum_d = |i: usize| -> f32 { active.iter().map(|&k| distances[i][k]).sum::<f32>() };

        // Find i,j for which Q(i,j) is minimal.
        let q = |&(&i, &j): &(&usize, &usize)| -> NotNan<f32> {
            let r = NotNan::new((active.len() - 2) as f32 * distances[i][j] - sum_d(i) - sum_d(j))
                .unwrap();
            r
        };

        // Find minimal distance pair.
        let (&i, &j) = active
            .iter()
            .cartesian_product(active.iter())
            .filter(|&(&i, &j)| i != j)
            .min_by_key(q)
            .unwrap();

        let (i, j) = (min(i, j), max(i, j));

        // Compute distance from merged vertex to the nodes being merged.
        let di = distances[i][j] / 2. + (sum_d(i) - sum_d(j)) / (2. * (active.len() as f32 - 2.));
        let dj = distances[i][j] - di;

        // Remove j from positions considered in later iterations.
        active.remove(active.iter().position(|&x| x == j).unwrap());

        // Compute all other distances.
        active.iter().filter(|&&k| k != i).for_each(|&k| {
            let dk = (distances[i][k] + distances[j][k] - distances[i][j]) / 2.;
            distances[i][k] = dk;
            distances[k][i] = dk;
        });
        parts[i] = Some(Phylogeny::join(
            parts[i].take().unwrap(),
            di,
            parts[j].take().unwrap(),
            dj,
        ));
    }

    // Merge the two remaining vertices with the given distance.
    if let [i, j] = active[..] {
        let d = distances[i][j] / 2.;
        parts[i] = Some(Phylogeny::join(
            parts[i].take().unwrap(),
            d,
            parts[j].take().unwrap(),
            d,
        ));
    }

    parts[0].take().unwrap()
}

/// Compute the Robingson-Foulds metric.
///
/// Also called the RF-distance, a distance metric between two phylogenetic trees.
///
/// Leaf nodes must be named. Names for internal nodes are ignored.
pub fn robinson_foulds_distance(p: &Phylogeny, q: &Phylogeny) -> usize {
    type Partition = BTreeSet<String>;
    type Partitions = HashSet<Partition>;

    fn leaf_dfs(node: &Phylogeny, leafs: &mut Partition) {
        if node.children.is_empty() {
            leafs.insert(node.name.to_string());
        } else {
            node.children.iter().for_each(|(c, _d)| leaf_dfs(c, leafs));
        }
    }
    fn all_leafs(p: &Phylogeny) -> Partition {
        let mut leafs = Partition::new();
        // Special case when the root is a leaf.
        if p.children.len() == 1 && !p.name.is_empty() {
            leafs.insert(p.name.to_string());
        }
        leaf_dfs(p, &mut leafs);
        leafs
    }
    let leafs_p = all_leafs(p);
    let leafs_q = all_leafs(q);
    assert_eq!(leafs_p, leafs_q);

    fn insert_part_and_complement(
        part: &BTreeSet<String>,
        complement: &BTreeSet<String>,
        parts: &mut HashSet<BTreeSet<String>>,
    ) {
        println!("insert {:?}", part);
        parts.insert(part.clone());
        let mut c = complement.clone();
        part.iter().for_each(|x| {
            c.remove(x);
        });
        parts.insert(c);
    }

    fn parts_dfs(node: &Phylogeny, complement: &Partition, parts: &mut Partitions) -> Partition {
        let mut part = Partition::new();
        if node.children.is_empty() {
            part.insert(node.name.to_string());
        } else {
            for (c, _d) in &node.children {
                for x in parts_dfs(&c, complement, parts) {
                    part.insert(x);
                }
            }
        }
        insert_part_and_complement(&part, complement, parts);
        return part;
    }
    let (mut partitions_p, mut partitions_q) = (Partitions::new(), Partitions::new());
    // Find the partitions of p and q.
    parts_dfs(p, &leafs_p, &mut partitions_p);
    parts_dfs(q, &leafs_p, &mut partitions_q);
    println!("Parts P: {:?}", partitions_p);
    println!("Parts Q: {:?}", partitions_q);
    if leafs_p.len() > 1 {
        assert_eq!(partitions_p.len(), 4 * leafs_p.len() - 4);
        assert_eq!(partitions_q.len(), 4 * leafs_p.len() - 4);
    }

    // Find the number of parts in one but not the other.
    let cnt = partitions_p
        .iter()
        .filter(|x| !partitions_q.contains(x))
        .count()
        + partitions_q
            .iter()
            .filter(|x| !partitions_p.contains(x))
            .count();
    println!("RF Distance {}", cnt);
    assert!(cnt % 2 == 0);
    cnt / 2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phylogeny_to_from_string() {
        // From Wikipedia: https://en.wikipedia.org/wiki/Newick_format
        let strings = vec![
            // This is not supported currently
            //"(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;",  //all have a distance to parent
            "(,,(,));",                            //no nodes are named
            "(A,B,(C,D));",                        //leaf nodes are named
            "(A,B,(C,D)E)F;",                      //all nodes are named
            "(:0.1,:0.2,(:0.3,:0.4):0.5);",        //all but root node have a distance to parent
            "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",    //distances and leaf names (popular)
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",  //distances and all names
            "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;", //a tree rooted on a leaf node (rare)
        ];
        for s in strings {
            assert_eq!(s.parse::<Phylogeny>().unwrap().to_string(), s);
        }
    }

    #[test]
    fn phylogeny_value() {
        let s = "(((a:8.5,b:8.5):2.5,e:11):5.5,(c:14,d:14):2.5);";
        let p = s.parse::<Phylogeny>().unwrap();
        assert_eq!(p.children[0].0.children[0].0.children[1].0.name, "b");
        assert_eq!(p.children[0].0.children[0].0.children[1].1, Some(8.5));
        assert_eq!(p.to_string(), s);
    }

    #[test]
    fn upgma_wiki() {
        // Sample from Wikipedia
        let d = DistanceMatrix::new(
            vec!["a", "b", "c", "d", "e"]
                .iter()
                .map(|&s| s.into())
                .collect(),
            vec![
                vec![0., 17., 21., 31., 23.],
                vec![17., 0., 30., 34., 21.],
                vec![32., 40., 0., 28., 39.],
                vec![31., 34., 28., 0., 43.],
                vec![23., 21., 39., 43., 0.],
            ],
        );
        let s = "(((a:8.5,b:8.5):2.5,e:11):5.5,(c:14,d:14):2.5);";
        assert_eq!(upgma(&d), s.parse().unwrap());
    }

    #[test]
    fn neighbor_joining_wiki() {
        // Sample from Wikipedia
        let d = DistanceMatrix::new(
            vec!["a", "b", "c", "d", "e"]
                .iter()
                .map(|&s| s.into())
                .collect(),
            vec![
                vec![0., 5., 9., 9., 8.],
                vec![5., 0., 10., 10., 9.],
                vec![9., 10., 0., 8., 7.],
                vec![9., 10., 8., 0., 3.],
                vec![8., 9., 7., 3., 0.],
            ],
        );
        let s = "((((a:2,b:3):3,c:4):2,d:2):0.5,e:0.5);";
        assert_eq!(neighbor_joining(&d), s.parse().unwrap());
    }

    #[test]
    fn robinson_foulds_distance_small() {
        let p1 = "a;".parse().unwrap();
        let p2 = "(a);".parse().unwrap();
        let p3 = "((a));".parse().unwrap();
        assert_eq!(robinson_foulds_distance(&p1, &p2), 0);
        assert_eq!(robinson_foulds_distance(&p1, &p3), 0);
        assert_eq!(robinson_foulds_distance(&p2, &p3), 0);

        let p1 = "(a,b,c);".parse().unwrap();
        let p2 = "((c,a),b);".parse().unwrap();
        let p3 = "(c,(b,a));".parse().unwrap();
        assert_eq!(robinson_foulds_distance(&p1, &p2), 0);
        assert_eq!(robinson_foulds_distance(&p1, &p3), 0);
        assert_eq!(robinson_foulds_distance(&p2, &p3), 0);

        let p1 = "((a,b),(c,d));".parse().unwrap();
        let p2 = "(((a,b),c),d);".parse().unwrap();
        let p3 = "((a,c),(b,d));".parse().unwrap();
        assert_eq!(robinson_foulds_distance(&p1, &p2), 0);
        assert_eq!(robinson_foulds_distance(&p1, &p3), 2);
        assert_eq!(robinson_foulds_distance(&p2, &p3), 2);
    }

    #[test]
    #[should_panic]
    fn robinson_foulds_distance_different_labels() {
        let p1 = "((a,b),(c,d));".parse().unwrap();
        let p2 = "(((a,B),c),d);".parse().unwrap();
        assert_eq!(robinson_foulds_distance(&p1, &p2), 0);
    }
}
