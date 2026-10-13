"""Driver scores and Pratt pair attribution, against the R package.

Both languages call the same C++, so a disagreement here is a marshalling
error -- wrong array order, wrong shape, a transposed matrix -- which is
exactly the class of bug worth a test.
"""
from pathlib import Path

import numpy as np

import pace

FIX = Path(__file__).parent / "fixtures"


def col(name):
    return np.ascontiguousarray(np.loadtxt(FIX / (name + ".gz")).ravel())


def test_driver_scores_match_r():
    out = pace.driver_scores(
        estimate_shrunk=col("ds_estimate.txt"),
        spec=col("ds_spec.txt"),
        focal_mean=col("ds_focal_mean.txt"),
        mu_bar=col("ds_mu_bar.txt"),
        alpha=col("ds_alpha.txt"),
        var_n=float(np.loadtxt(FIX / "ds_var_n.txt.gz")),
    )
    # the R binding returns the core's own names; ours are the renamed ones
    for key, fixture in [
        ("driver_score", "ds_expect_MCSD.txt"),
        ("driver_score_spec4", "ds_expect_MCSD4.txt"),
        ("V_S", "ds_expect_V_S.txt"),
        ("V_total", "ds_expect_V_total.txt"),
        ("R2_S", "ds_expect_R2_S.txt"),
        ("V_resid", "ds_expect_V_resid.txt"),
    ]:
        np.testing.assert_allclose(out[key], col(fixture), rtol=1e-10, atol=0,
                                   err_msg=f"{key} disagrees with the R package")


def test_pratt_matches_r():
    n_pairs = 7
    sigma = np.asfortranarray(np.loadtxt(FIX / "pratt_sigma.txt.gz").reshape(n_pairs, n_pairs))
    u = np.asfortranarray(np.loadtxt(FIX / "pratt_u.txt.gz").reshape(n_pairs, -1))
    out = pace.pair_variance_pratt(sigma, u)
    for key, fixture in [
        ("v_pair", "pratt_expect_V_pair_pratt.txt"),
        ("v_diag", "pratt_expect_V_pair_diag.txt"),
        ("share_pct", "pratt_expect_within_focal_share_pct.txt"),
        ("diag_pct", "pratt_expect_within_focal_diag_pct.txt"),
        ("sign", "pratt_expect_sign.txt"),
    ]:
        np.testing.assert_allclose(out[key], col(fixture), rtol=1e-10, atol=1e-300,
                                   err_msg=f"{key} disagrees with the R package")
    # the scalars R reports alongside
    assert out["n_negative"] == int(np.loadtxt(FIX / "pratt_expect_n_negative_pairs.txt.gz"))
    np.testing.assert_allclose(out["cross_cov_pct"],
                               float(np.loadtxt(FIX / "pratt_expect_cross_cov_pct.txt.gz")),
                               rtol=1e-10)
    # shares are a within-focal attribution, so they sum to 100
    np.testing.assert_allclose(out["share_pct"].sum(), 100.0, rtol=1e-8)
