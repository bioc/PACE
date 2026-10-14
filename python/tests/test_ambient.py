"""The ambient contamination field, against the R package.

Compared on real coordinates: the sparsity pattern as well as the values,
since a field with the right numbers in the wrong places would still look
reasonable.
"""
from pathlib import Path

import numpy as np

import pace

FIX = Path(__file__).parent / "fixtures"


def col(name):
    return np.loadtxt(FIX / (name + ".gz"))


def build():
    x = np.ascontiguousarray(col("kern_x.txt"))
    y = np.ascontiguousarray(col("kern_y.txt"))
    cell_type = col("kern_type.txt").astype(np.int32)
    image = np.zeros(x.size, dtype=np.int32)
    return pace.ambient_field(x, y, cell_type, image, n_images=1, h_tech=5.0,
                              edge_correct=True)


def test_ambient_field_matches_r():
    out = build()
    # the sparsity pattern, not just the values
    np.testing.assert_array_equal(out["indices"], col("amb_expect_i.txt").astype(np.int32))
    np.testing.assert_array_equal(out["indptr"], col("amb_expect_p.txt").astype(np.int64))
    np.testing.assert_allclose(out["data"], col("amb_expect_x.txt"), rtol=1e-12, atol=0)
    np.testing.assert_allclose(out["edge_fraction"], col("amb_expect_edge.txt"),
                               rtol=1e-12, atol=0)


def test_the_field_has_the_properties_it_should():
    out = build()
    n = out["shape"][0]
    assert out["indptr"].size == n + 1
    assert out["indptr"][0] == 0 and out["indptr"][-1] == out["nnz"]
    assert out["data"].size == out["nnz"] == out["indices"].size
    # weights are exp(-d/h) over an edge fraction, so strictly positive
    assert out["data"].min() > 0
    # edge fractions are a fraction of a disc, and 1 where no correction applies
    assert out["edge_fraction"].min() > 0
    assert out["edge_fraction"].max() <= 1.0
    # no cell contaminates itself
    for j in range(0, n, 500):
        rows = out["indices"][out["indptr"][j]:out["indptr"][j + 1]]
        assert j not in rows


def test_turning_off_the_edge_correction_sets_every_fraction_to_one():
    x = np.ascontiguousarray(col("kern_x.txt"))
    y = np.ascontiguousarray(col("kern_y.txt"))
    cell_type = col("kern_type.txt").astype(np.int32)
    image = np.zeros(x.size, dtype=np.int32)
    out = pace.ambient_field(x, y, cell_type, image, n_images=1, h_tech=5.0,
                             edge_correct=False)
    np.testing.assert_allclose(out["edge_fraction"], 1.0)
