"""Builds the Python bindings against the SAME C++ core the R package uses.

Nothing here is a reimplementation: ../src/core/mash.cpp is compiled directly,
so the two languages cannot drift apart. The core includes no R headers, which
is what makes this possible.
"""
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# setuptools rejects absolute paths in `sources`, so the core is referenced
# relative to this file. It still compiles the one shared copy; nothing is
# duplicated into the Python tree.
CORE = Path("..") / "src" / "core"

ext_modules = [
    Pybind11Extension(
        "pace._core",
        sources=["src/pace_py.cpp", str(CORE / "mash.cpp")],
        include_dirs=[str(CORE)],
        cxx_std=17,
    )
]

setup(cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)
