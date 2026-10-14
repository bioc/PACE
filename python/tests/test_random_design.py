"""The random-effect design block, against the R package.

Checks the sparse structure and the grouping, not just the values: a block
with the right numbers in the wrong columns would fit a different model.
"""
from pathlib import Path

import numpy as np

import pace

FIX = Path(__file__).parent / "fixtures"


def col(name):
    return np.loadtxt(FIX / (name + ".gz"))


def inputs():
    n, n_terms, n_groups = (int(v) for v in col("rd_dims.txt"))
    x_terms = np.asfortranarray(col("rd_x.txt").reshape(n, n_terms, order="F"))
    cell_group = col("rd_group.txt").astype(np.int32)
    return x_terms, cell_group, n, n_terms, n_groups


def test_design_block_matches_r():
    x_terms, cell_group, n, n_terms, n_groups = inputs()
    out = pace.random_design_block(x_terms, cell_group, n_groups)

    assert out["shape"] == (n, n_terms * n_groups)
    np.testing.assert_array_equal(out["indices"], col("rd_expect_i.txt").astype(np.int32))
    np.testing.assert_array_equal(out["indptr"], col("rd_expect_p.txt").astype(np.int32))
    np.testing.assert_allclose(out["data"], col("rd_expect_x.txt"), rtol=1e-14, atol=0)


def test_the_grouping_matches_r():
    x_terms, cell_group, n, n_terms, n_groups = inputs()
    out = pace.random_design_block(x_terms, cell_group, n_groups)
    np.testing.assert_array_equal(out["group_start"], col("rd_expect_gstart.txt").astype(np.int32))
    np.testing.assert_array_equal(out["cells_by_group"], col("rd_expect_cbg.txt").astype(np.int32))


def test_the_grouping_partitions_the_cells():
    x_terms, cell_group, n, n_terms, n_groups = inputs()
    out = pace.random_design_block(x_terms, cell_group, n_groups)
    start, cbg = out["group_start"], out["cells_by_group"]

    assert start[0] == 0 and start[-1] == n
    assert sorted(cbg.tolist()) == list(range(n))   # every cell once
    for g in range(n_groups):
        members = cbg[start[g]:start[g + 1]]
        assert np.all(cell_group[members] == g)     # and in the right group
        assert np.all(np.diff(members) > 0)         # ascending, as split() gives
