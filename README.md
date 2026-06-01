# rsomics-fcluster

Form flat clusters from a hierarchical linkage matrix by cutting it — a Rust
reimplementation of `scipy.cluster.hierarchy.fcluster`. Reads the linkage
matrix written by [`rsomics-upgma`](https://github.com/omics-rust/rsomics-upgma)
(or scipy's `linkage`) and emits one integer cluster label per observation,
byte-for-byte matching scipy.

## Usage

```
rsomics-fcluster <linkage.tsv> --threshold T [--criterion C] [-o labels.tsv]
```

The linkage matrix is `n-1` tab-separated rows of `left right height size`;
cluster ids `< n` are original observations, ids `>= n` reference the cluster
formed at row `id - n`. Output is `n` lines, label for observation `0..n`.

| Criterion | Meaning |
|---|---|
| `distance` (default) | Cut so no two observations in a cluster have cophenetic distance `> T`. |
| `maxclust` | Smallest threshold yielding at most `T` clusters. |
| `inconsistent` | Cut where the inconsistency coefficient (over `--depth` levels, default 2) is `<= T`. |
| `monocrit` | Cut where a user-supplied per-node statistic (`--monocrit FILE`, `n-1` values) is `<= T`. |

```
rsomics-fcluster Z.tsv --threshold 1.1 --criterion distance
rsomics-fcluster Z.tsv --threshold 4   --criterion maxclust
rsomics-fcluster Z.tsv --threshold 0.9 --criterion inconsistent --depth 3
```

## Origin

This crate is an independent Rust reimplementation of
`scipy.cluster.hierarchy.fcluster`, based on reading and citing SciPy's
`_hierarchy.pyx` (`cluster_dist`, `cluster_maxclust_dist`, `cluster_in`,
`cluster_monocrit`, `cluster_maxclust_monocrit`, `get_max_dist_for_each_cluster`,
`get_max_Rfield_for_each_cluster`, `inconsistent`) — BSD-3-Clause.

The cluster-label numbering reproduces scipy's exact DFS traversal in
`cluster_monocrit` (push the left child first; clusters are numbered in the
order leaves are reached on the descent), so output integers match scipy
element-for-element. The maxclust binary search mirrors
`cluster_maxclust_monocrit` (search the sorted cophenetic thresholds for the
smallest one giving `<= t` clusters).

License: MIT OR Apache-2.0.
Upstream credit: SciPy (https://scipy.org, BSD-3-Clause).
