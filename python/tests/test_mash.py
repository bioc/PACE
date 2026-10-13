"""The Python bindings must reproduce the R package's published call counts.

The fixtures are every neighbour slice of both canonical cohorts, dumped with
mashr's own lfsr alongside. The target is exact: 1634 breast cancer calls and
426 melanoma calls, which is what the manuscript reports.
"""
from pathlib import Path

import numpy as np
import pytest

import pace

SLICES = Path(__file__).parent / "fixtures"
EXPECTED = {"bc": 1634, "mel": 426}


def load(stem):
    """mashr wrote these row-major; numpy reads them that way too."""
    n_conditions = int(np.loadtxt(SLICES / f"{stem}_dim.txt.gz"))
    bhat = np.loadtxt(SLICES / f"{stem}_b.txt.gz").reshape(-1, n_conditions)
    shat = np.loadtxt(SLICES / f"{stem}_s.txt.gz").reshape(-1, n_conditions)
    lfsr = np.loadtxt(SLICES / f"{stem}_lfsr.txt.gz").reshape(-1, n_conditions)
    return np.asfortranarray(bhat), np.asfortranarray(shat), lfsr


def slices_for(cohort):
    return sorted(p.name[: -len("_b.txt.gz")] for p in SLICES.glob(f"{cohort}_*_b.txt.gz"))


@pytest.mark.parametrize("cohort", ["bc", "mel"])
def test_reproduces_published_call_counts(cohort):
    total = 0
    for stem in slices_for(cohort):
        bhat, shat, lfsr_r = load(stem)
        out = pace.mash(bhat, shat)
        called = out["lfsr"] < 0.05
        reference = lfsr_r < 0.05
        # not merely the same count: the same genes
        assert np.array_equal(called, reference), f"{stem} disagrees with mashr"
        total += int(called.sum())
    assert total == EXPECTED[cohort]


def test_posterior_shapes_and_finiteness():
    bhat, shat, _ = load("bc_Dendritic_Cell")
    out = pace.mash(bhat, shat)
    for key in ("posterior_mean", "posterior_sd", "lfsr"):
        assert out[key].shape == bhat.shape
        assert np.all(np.isfinite(out[key]))
    assert np.all(out["posterior_sd"] >= 0)
    assert np.all((out["lfsr"] >= 0) & (out["lfsr"] <= 1))
    assert abs(sum(out["pi"]) - 1.0) < 1e-10
