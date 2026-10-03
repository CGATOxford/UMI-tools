"""Regression tests for neighbour enumeration and adjacency clustering."""

import itertools
import random

import pytest

from umi_tools import network


def original_neighbours(umis, substr_idx):
    """Reference implementation before incremental prefix tracking."""
    for i, umi in enumerate(umis, 1):
        neighbours = set()
        for idx, substr_map in substr_idx.items():
            neighbours = neighbours.union(substr_map[umi[slice(*idx)]])
        neighbours.difference_update(umis[:i])
        for neighbour in neighbours:
            yield umi, neighbour


def original_min_account(self, cluster, adj_list, counts):
    if len(cluster) == 1:
        return list(cluster)
    nodes = sorted(cluster, key=lambda x: counts[x], reverse=True)
    for i in range(len(nodes) - 1):
        if not network.remove_umis(adj_list, cluster, nodes[:i+1]):
            return nodes[:i+1]


def make_counts(seed, size, length, alphabet=b"ACGT"):
    rng = random.Random(seed)
    counts = {}
    while len(counts) < size:
        umi = bytes(rng.choices(alphabet, k=length))
        counts[umi] = rng.choice([1, 1, 2, 3, 10, 20])
    return counts


@pytest.mark.parametrize("seed", range(8))
@pytest.mark.parametrize("length,size,threshold,alphabet", [
    (12, 1000, 1, b"ACGT"),
    (8, 200, 2, b"ACGT"),
    (4, 200, 1, b"ACGTN"),
    (3, 64, 3, b"ACGT"),
    (1, 4, 2, b"ACGT"),
    (6, 100, 0, b"ACGT"),
])
def test_neighbour_pairs(seed, length, size, threshold, alphabet):
    umis = list(make_counts(seed, size, length, alphabet))
    index = network.build_substr_idx(umis, length, threshold)
    actual = list(network.iter_nearest_neighbours(umis, index))
    assert actual == list(original_neighbours(umis, index))
    assert len(actual) == len(set(actual))
    close_pairs = {pair for pair in itertools.combinations(umis, 2)
                   if network.edit_distance(*pair) <= threshold}
    assert close_pairs <= set(actual)


@pytest.mark.parametrize("method", ["directional", "adjacency", "cluster"])
@pytest.mark.parametrize("seed", range(8))
@pytest.mark.parametrize("length,size,threshold", [
    (8, 1, 1), (8, 25, 1), (8, 26, 1), (12, 1000, 1),
    (4, 200, 1), (6, 200, 2), (3, 64, 3), (6, 100, 0),
])
def test_cluster_output(monkeypatch, method, seed, length, size, threshold):
    counts = make_counts(seed, size, length)
    actual = network.UMIClusterer(method)(counts, threshold)
    with monkeypatch.context() as patch:
        patch.setattr(network, "iter_nearest_neighbours", original_neighbours)
        patch.setattr(network.UMIClusterer, "_get_best_min_account",
                      original_min_account)
        expected = network.UMIClusterer(method)(counts, threshold)
    # Preserve representatives, ties, group order and member order.
    assert actual == expected


@pytest.mark.parametrize("size", [1, 2, 3, 50, 500])
@pytest.mark.parametrize("shape", ["path", "star", "complete", "disconnected"])
def test_adjacency_cover(size, shape):
    nodes = list(range(size))
    if shape == "path":
        graph = {i: [j for j in (i-1, i+1) if 0 <= j < size]
                 for i in nodes}
    elif shape == "star":
        graph = {i: ([0] if i else nodes[1:]) for i in nodes}
    elif shape == "complete":
        graph = {i: [j for j in nodes if j != i] for i in nodes}
    else:
        graph = {i: [] for i in nodes}
    for counts in ({i: 1 for i in nodes}, {i: i+1 for i in nodes}):
        clusterer = network.UMIClusterer("adjacency")
        assert clusterer._get_best_min_account(nodes, graph, counts) == \
            original_min_account(clusterer, nodes, graph, counts)
