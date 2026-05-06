# Release Structure

This repository is organized by scientific maturity rather than by the order in
which the material was created.

## 2021 Single-Phase Material

`docs/2021-single-phase/` preserves the explanatory notes and examples
associated with the original single-phase theoretical publication.  This folder
is intentionally historical: it records the starting point of the method and is
not the active implementation used for the current three-phase studies.

## Current MATLAB API

`matlab/+spacor/` contains the public MATLAB namespace used by the current
reproducibility package.  It exposes stable entry points for:

- synthetic fixture generation,
- single-phase load-equivalent characterization,
- three-wire load-equivalent characterization,
- finite-window campaigns,
- measurement-degradation campaigns,
- derivative and primitive reconstruction.

The implementation emphasizes reproducible sampled-data identification:
finite windows, matrix estimators, harmonic projection for derivatives and
integrals, information diagnostics, residual checks, and passive model
selection.

## Validated Runtime

`matlab/legacy/runtime/` contains validated research routines that are still
used internally by the public API.  They are kept separate because they preserve
the numerical behavior used in the paper-development campaigns, while the
`+spacor` namespace provides the stable user-facing layer.

## Public Experiments

`scripts/` and `tests/` provide entry points for regenerating synthetic results.
They do not require measurement datasets.  Generated outputs are written under
ignored result folders and are not committed to Git.

## Measurement Data

No private measurement data are shipped in this release.  Any future public
dataset should be released separately with clear consent, metadata, license,
and integrity checks.
