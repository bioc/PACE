"""PACE: Proximity-Associated Changes in Expression.

This package binds the same compiled core the R package uses, so results do
not depend on which language you call it from.
"""

from ._core import (
    ambient_field,
    driver_scores,
    mash,
    neighbour_kernels,
    pair_variance_pratt,
)

__all__ = [
    "mash",
    "driver_scores",
    "pair_variance_pratt",
    "neighbour_kernels",
    "ambient_field",
]
__version__ = "0.0.1"
