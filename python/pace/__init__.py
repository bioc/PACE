"""PACE: Proximity-Associated Changes in Expression.

This package binds the same compiled core the R package uses, so results do
not depend on which language you call it from.
"""

from ._core import mash

__all__ = ["mash"]
__version__ = "0.0.1"
