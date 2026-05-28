# TODO: Install prolfquabenchmark without vignettes

## Problem

`prolfquabenchmark` vignettes are the main content and can be expensive to build. For local ecosystem install chains we
want `make install` for this package only to install without building vignettes, while leaving the shared R package
Makefile template and other package install targets unchanged.

## Plan

- Keep `make build`, `make build-vignettes`, and `make check` unchanged so explicit vignette builds/checks remain
  available.
- Add a package-local no-vignettes build command in `prolfquabenchmark/Makefile`.
- Change only the `install` target in `prolfquabenchmark/Makefile` to document, build without vignettes, then install the
  generated tarball.
- Update the local help text so the package-specific behavior is visible.
- Verify with `make -n install` that no shared template or sibling package behavior changed.
