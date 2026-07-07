use crate::linkage::Linkage;

/// Flat-cluster formation criterion. scipy `_hierarchy.pyx`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Criterion {
    Distance,
    MaxClust,
    Inconsistent,
    Monocrit,
}

/// Cophenetic max distance below and including each non-singleton node.
fn max_dist_for_each_cluster(z: &Linkage) -> Vec<f64> {
    monocrit_postorder(z, |root| z.height[root])
}

/// Max inconsistency coefficient (R column 3) below and including each node.
fn max_rfield_for_each_cluster(z: &Linkage, r3: &[f64]) -> Vec<f64> {
    monocrit_postorder(z, |root| r3[root])
}

/// Post-order over the tree taking, for each node, the max of its own statistic
/// and its children's. `self_stat(root)` gives the node's own value.
fn monocrit_postorder<F: Fn(usize) -> f64>(z: &Linkage, self_stat: F) -> Vec<f64> {
    let n = z.n;
    let mut md = vec![0.0f64; n - 1];
    let mut visited = vec![false; 2 * n - 1];
    let mut stack = vec![2 * n - 2];

    while let Some(&top) = stack.last() {
        let root = top - n;
        let i_lc = z.left[root];
        let i_rc = z.right[root];

        if i_lc >= n && !visited[i_lc] {
            visited[i_lc] = true;
            stack.push(i_lc);
            continue;
        }
        if i_rc >= n && !visited[i_rc] {
            visited[i_rc] = true;
            stack.push(i_rc);
            continue;
        }

        let mut m = self_stat(root);
        if i_lc >= n {
            m = m.max(md[i_lc - n]);
        }
        if i_rc >= n {
            m = m.max(md[i_rc - n]);
        }
        md[root] = m;
        stack.pop();
    }
    md
}

/// scipy `cluster_monocrit`: cut where the monotone criterion `mc[node] <=
/// cutoff` first holds on a root-to-leaf descent, numbering clusters in the DFS
/// order leaves are reached. Produces scipy's exact integer labels.
fn cluster_monocrit(z: &Linkage, mc: &[f64], cutoff: f64) -> Vec<i32> {
    let n = z.n;
    let mut t = vec![0i32; n];
    let mut visited = vec![false; 2 * n - 1];
    let mut curr = vec![0usize; n];

    let mut n_cluster = 0i32;
    let mut cluster_leader: i64 = -1;
    let mut k: i64 = 0;
    curr[0] = 2 * n - 2;

    while k >= 0 {
        let root = curr[k as usize] - n;
        let i_lc = z.left[root];
        let i_rc = z.right[root];

        if cluster_leader == -1 && mc[root] <= cutoff {
            cluster_leader = root as i64;
            n_cluster += 1;
        }

        if i_lc >= n && !visited[i_lc] {
            visited[i_lc] = true;
            k += 1;
            curr[k as usize] = i_lc;
            continue;
        }
        if i_rc >= n && !visited[i_rc] {
            visited[i_rc] = true;
            k += 1;
            curr[k as usize] = i_rc;
            continue;
        }

        if i_lc < n {
            if cluster_leader == -1 {
                n_cluster += 1;
            }
            t[i_lc] = n_cluster;
        }
        if i_rc < n {
            if cluster_leader == -1 {
                n_cluster += 1;
            }
            t[i_rc] = n_cluster;
        }

        if cluster_leader == root as i64 {
            cluster_leader = -1;
        }
        k -= 1;
    }
    t
}

/// scipy `cluster_maxclust_monocrit`: binary-search the smallest monocrit
/// threshold yielding `<= max_nc` clusters, then label via `cluster_monocrit`.
fn cluster_maxclust_monocrit(z: &Linkage, mc: &[f64], max_nc: usize) -> Vec<i32> {
    let n = z.n;
    if max_nc >= n {
        return (1..=n as i32).collect();
    }
    // scipy's binary search leaves upper_idx at n-1 when no split can honour a
    // <=0 cluster budget, then reads MC[n-1] one past the length-(n-1) array —
    // undefined behaviour that happens to yield the all-singletons partition.
    // fcluster(Z, 0, 'maxclust') is documented to place every point in its own
    // cluster, so form that DFS-ordered labelling directly (a cutoff below every
    // monocrit value makes each leaf its own cluster).
    if max_nc == 0 {
        return cluster_monocrit(z, mc, f64::NEG_INFINITY);
    }

    let mut visited = vec![false; 2 * n - 1];
    let mut curr = vec![0usize; n];

    let mut lower_idx: i64 = -1;
    let mut upper_idx: usize = n - 1;

    while upper_idx as i64 - lower_idx > 1 {
        let i = ((lower_idx + upper_idx as i64) >> 1) as usize;
        let thresh = mc[i];

        visited.iter_mut().for_each(|v| *v = false);
        let mut nc = 0usize;
        let mut k: i64 = 0;
        curr[0] = 2 * n - 2;

        while k >= 0 {
            let root = curr[k as usize] - n;
            let i_lc = z.left[root];
            let i_rc = z.right[root];

            if mc[root] <= thresh {
                nc += 1;
                if nc > max_nc {
                    break;
                }
                k -= 1;
                visited[i_lc] = true;
                visited[i_rc] = true;
                continue;
            }

            if !visited[i_lc] {
                visited[i_lc] = true;
                if i_lc >= n {
                    k += 1;
                    curr[k as usize] = i_lc;
                    continue;
                }
                nc += 1;
                if nc > max_nc {
                    break;
                }
            }

            if !visited[i_rc] {
                visited[i_rc] = true;
                if i_rc >= n {
                    k += 1;
                    curr[k as usize] = i_rc;
                    continue;
                }
                nc += 1;
                if nc > max_nc {
                    break;
                }
            }

            k -= 1;
        }

        if nc > max_nc {
            lower_idx = i as i64;
        } else {
            upper_idx = i;
        }
    }

    cluster_monocrit(z, mc, mc[upper_idx])
}

/// scipy `inconsistent`: for each non-singleton, mean/std of merge heights over
/// the link plus its descendants down to `depth` levels; column 3 is the
/// inconsistency coefficient `(Z[i].height - mean) / std`.
fn inconsistent_r3(z: &Linkage, depth: usize) -> Vec<f64> {
    let n = z.n;
    let d = depth.max(1);
    let mut r3 = vec![0.0f64; n - 1];
    let mut visited = vec![false; 2 * n - 1];
    let mut curr = vec![0usize; n];

    // `i` is the node index — it seeds the traversal and indexes z.height/r3.
    #[allow(clippy::needless_range_loop)]
    for i in 0..n - 1 {
        let mut k: i64 = 0;
        let mut level_count = 0u64;
        let mut level_sum = 0.0f64;
        let mut level_std_sum = 0.0f64;
        visited.iter_mut().for_each(|v| *v = false);
        curr[0] = i;

        while k >= 0 {
            let root = curr[k as usize];

            if (k as usize) < d - 1 {
                let i_lc = z.left[root];
                if i_lc >= n && !visited[i_lc] {
                    visited[i_lc] = true;
                    k += 1;
                    curr[k as usize] = i_lc - n;
                    continue;
                }
                let i_rc = z.right[root];
                if i_rc >= n && !visited[i_rc] {
                    visited[i_rc] = true;
                    k += 1;
                    curr[k as usize] = i_rc - n;
                    continue;
                }
            }

            let dist = z.height[root];
            level_count += 1;
            level_sum += dist;
            level_std_sum += dist * dist;
            k -= 1;
        }

        let lc = level_count as f64;
        let mean = level_sum / lc;
        let level_std = if level_count < 2 {
            (level_std_sum - level_sum * level_sum) / lc
        } else {
            (level_std_sum - (level_sum * level_sum) / lc) / (lc - 1.0)
        };
        if level_std > 0.0 {
            r3[i] = (z.height[i] - mean) / level_std.sqrt();
        }
    }
    r3
}

/// Form flat clusters. `t` is a cophenetic/inconsistency threshold for
/// distance/inconsistent/monocrit, or the max cluster count for maxclust.
/// `depth` only applies to inconsistent; `monocrit` only to the monocrit
/// criterion (the per-node statistic vector, length `n-1`).
#[must_use]
pub fn fcluster(
    z: &Linkage,
    t: f64,
    criterion: Criterion,
    depth: usize,
    monocrit: Option<&[f64]>,
) -> Vec<i32> {
    match criterion {
        Criterion::Distance => {
            let md = max_dist_for_each_cluster(z);
            cluster_monocrit(z, &md, t)
        }
        Criterion::Inconsistent => {
            let r3 = inconsistent_r3(z, depth);
            let mr = max_rfield_for_each_cluster(z, &r3);
            cluster_monocrit(z, &mr, t)
        }
        Criterion::Monocrit => cluster_monocrit(z, monocrit.unwrap(), t),
        Criterion::MaxClust => {
            let md = max_dist_for_each_cluster(z);
            cluster_maxclust_monocrit(z, &md, t as usize)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // scipy ward(pdist(X)) for the 12-point doc example; heights collapsed to
    // their exact computed values.
    fn ward_example() -> Linkage {
        let h12 = 1.2071067811865475;
        let h6 = 3.396751841199101;
        let h_root = 4.092065228194317;
        Linkage {
            left: vec![0, 3, 6, 9, 2, 5, 8, 11, 16, 18, 20],
            right: vec![1, 4, 7, 10, 12, 13, 14, 15, 17, 19, 21],
            height: vec![1.0, 1.0, 1.0, 1.0, h12, h12, h12, h12, h6, h6, h_root],
            n: 12,
        }
    }

    #[test]
    fn distance_matches_scipy() {
        let z = ward_example();
        assert_eq!(
            fcluster(&z, 2.0, Criterion::Distance, 2, None),
            vec![1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
        );
        assert_eq!(
            fcluster(&z, 4.0, Criterion::Distance, 2, None),
            vec![1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
        );
        // t larger than the root height: one cluster.
        assert_eq!(fcluster(&z, 9.0, Criterion::Distance, 2, None), vec![1; 12]);
        // t below every merge: every point its own cluster.
        assert_eq!(
            fcluster(&z, 0.5, Criterion::Distance, 2, None),
            (1..=12).collect::<Vec<_>>()
        );
    }

    #[test]
    fn maxclust_matches_scipy() {
        let z = ward_example();
        assert_eq!(
            fcluster(&z, 4.0, Criterion::MaxClust, 2, None),
            vec![1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
        );
        assert_eq!(
            fcluster(&z, 3.0, Criterion::MaxClust, 2, None),
            vec![1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
        );
        // max_nc >= n short-circuits to singletons.
        assert_eq!(
            fcluster(&z, 12.0, Criterion::MaxClust, 2, None),
            (1..=12).collect::<Vec<_>>()
        );
    }

    #[test]
    fn maxclust_zero_is_all_singletons() {
        let z = ward_example();
        // scipy fcluster(Z, 0, 'maxclust') -> every point its own cluster,
        // numbered in DFS order (leaf order here). Must not panic.
        assert_eq!(
            fcluster(&z, 0.0, Criterion::MaxClust, 2, None),
            (1..=12).collect::<Vec<_>>()
        );
        // A negative threshold saturates to max_nc = 0 through `t as usize`;
        // scipy likewise yields the degenerate all-singletons partition.
        assert_eq!(
            fcluster(&z, -1.0, Criterion::MaxClust, 2, None),
            (1..=12).collect::<Vec<_>>()
        );
    }

    #[test]
    fn inconsistent_matches_scipy() {
        let z = ward_example();
        assert_eq!(
            fcluster(&z, 1.0, Criterion::Inconsistent, 2, None),
            vec![1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
        );
        assert_eq!(
            fcluster(&z, 0.7, Criterion::Inconsistent, 3, None),
            vec![1, 1, 2, 3, 3, 4, 5, 5, 6, 7, 7, 8]
        );
    }
}
