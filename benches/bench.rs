use std::hint::black_box;

use criterion::{Criterion as Bench, criterion_group, criterion_main};
use rsomics_fcluster::{Criterion, Linkage, fcluster};

/// A perfectly-balanced binary linkage over `n = 2^levels` leaves, heights
/// increasing by level — enough structure to exercise every traversal.
fn balanced_linkage(levels: u32) -> Linkage {
    let n = 1usize << levels;
    let mut left = Vec::with_capacity(n - 1);
    let mut right = Vec::with_capacity(n - 1);
    let mut height = Vec::with_capacity(n - 1);

    let mut nodes: Vec<usize> = (0..n).collect();
    let mut next = n;
    let mut h = 1.0f64;
    while nodes.len() > 1 {
        let mut up = Vec::with_capacity(nodes.len() / 2);
        for pair in nodes.chunks(2) {
            left.push(pair[0]);
            right.push(pair[1]);
            height.push(h);
            up.push(next);
            next += 1;
        }
        nodes = up;
        h += 1.0;
    }
    Linkage {
        left,
        right,
        height,
        n,
    }
}

fn bench(c: &mut Bench) {
    let z = balanced_linkage(13); // 8192 leaves
    c.bench_function("distance", |b| {
        b.iter(|| black_box(fcluster(&z, 5.0, Criterion::Distance, 2, None)))
    });
    c.bench_function("maxclust", |b| {
        b.iter(|| black_box(fcluster(&z, 64.0, Criterion::MaxClust, 2, None)))
    });
    c.bench_function("inconsistent", |b| {
        b.iter(|| black_box(fcluster(&z, 1.0, Criterion::Inconsistent, 2, None)))
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
