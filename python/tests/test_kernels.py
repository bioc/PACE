"""Neighbour kernels, against the R package on real coordinates.

A kernel built on a transposed coordinate matrix, or on 1-based type codes,
would be wrong in a way nothing downstream announces -- so this compares
against R on the packaged breast cancer subset rather than on a toy.
"""
from pathlib import Path

import numpy as np

import pace

FIX = Path(__file__).parent / "fixtures"


def col(name):
    return np.loadtxt(FIX / (name + ".gz"))


def test_kernels_match_r():
    x = np.ascontiguousarray(col("kern_x.txt"))
    y = np.ascontiguousarray(col("kern_y.txt"))
    cell_type = col("kern_type.txt").astype(np.int32)
    n_types = int(col("kern_ntypes.txt"))
    group = np.zeros(x.size, dtype=np.int32)

    out = pace.neighbour_kernels(x, y, cell_type, n_types, group,
                                 per_group=False, h_bio=30.0, h_tech=5.0, eps=90.0)

    for key, fixture in [("k_bio", "kern_expect_bio.txt"), ("k_tech", "kern_expect_tech.txt")]:
        expected = col(fixture).reshape(x.size, n_types)
        np.testing.assert_allclose(out[key], expected, rtol=1e-12, atol=1e-12,
                                   err_msg=f"{key} disagrees with the R package")


def test_kernels_have_the_shape_and_sign_they_should():
    x = np.ascontiguousarray(col("kern_x.txt"))
    y = np.ascontiguousarray(col("kern_y.txt"))
    cell_type = col("kern_type.txt").astype(np.int32)
    n_types = int(col("kern_ntypes.txt"))
    group = np.zeros(x.size, dtype=np.int32)
    out = pace.neighbour_kernels(x, y, cell_type, n_types, group)

    assert out["k_bio"].shape == (x.size, n_types)
    assert out["k_tech"].shape == (x.size, n_types)
    # kernel weights are counts of neighbours, so never negative
    assert out["k_bio"].min() >= 0
    assert out["k_tech"].min() >= 0
    # the default radius is three biological bandwidths
    assert out["eps"] == 90.0
    # the technical kernel is the tighter one, so it sees less
    assert out["k_tech"].sum() < out["k_bio"].sum()
