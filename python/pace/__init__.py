"""PACE: Proximity-Associated Changes in Expression.

This package binds the same compiled core the R package uses, so results do
not depend on which language you call it from.
"""

from ._core import driver_scores, mash, neighbour_kernels, pair_variance_pratt

__all__ = ["mash", "driver_scores", "pair_variance_pratt", "neighbour_kernels"]
__version__ = "0.0.1"
