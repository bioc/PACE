# PACE for Python

Bindings for the PACE compiled core. Nothing here reimplements the method: the
extension compiles `../src/core/*.cpp` directly, the same translation units the
R package uses, so the two cannot drift apart.

## Install

```sh
pip install -e python/
```

### macOS: pin the toolchain

A conda or Homebrew `clang` usually brings an old `ld64` along with it, and an
old `ld64` links a shared object that loads nowhere:

```
ImportError: dlopen(...): mis-aligned LINKEDIT string pool
```

The build succeeds, so this only shows up at import. Use Apple's toolchain:

```sh
CC=/usr/bin/clang CXX=/usr/bin/clang++ \
  LDSHARED="/usr/bin/clang -bundle -undefined dynamic_lookup" \
  pip install -e python/
```

`ld -v` should report `PROJECT:ld-27037` or newer. `ld64-530` is the one that
fails. The R package has the same requirement for the same reason.

## What is bound so far

| function | does |
|---|---|
| `mash` | canonical multivariate adaptive shrinkage for one neighbour slice |
| `driver_scores` | per-gene driver scores for one focal-neighbour pair |
| `pair_variance_pratt` | Pratt attribution of spatial variance across neighbours |
| `neighbour_kernels` | Gaussian and technical neighbour kernels from coordinates |

Fitting is still R only. `neighbour_kernels()` is its first stage and runs
standalone, but the rest -- the design matrices, the ambient field and the
IRLS loop -- is not bound yet, so you cannot go from an `AnnData` to a fit.
The other three functions take fitted output and summarise it.

## Tests

```sh
python -m pytest python/tests
```

The mash tests hold the bindings to the same target as the R package: every
neighbour slice of both canonical cohorts, compared gene by gene against mashr
rather than by count, giving 1634 breast cancer and 426 melanoma calls.
