# SPACOR

Geometric time-domain load-equivalent identification from terminal voltage and
current measurements.

This repository is being organized as the public reproducibility package for
the SPACOR line of work.  It separates the original 2021 single-phase material
from the newer MATLAB framework used for three-phase and sampled-data
identification studies.

## Repository Layout

```text
docs/2021-single-phase/   Original 2021 single-phase explanatory notes.
matlab/+spacor/           Current public MATLAB API.
matlab/legacy/runtime/    Validated research runtime used by the API.
scripts/                  Public smoke and reproduction entry points.
tests/                    Synthetic regression tests.
docs/                     Public API, campaigns, and data policy.
```

The files under `docs/2021-single-phase/` are preserved because they document
the original theory and examples associated with the 2021 TPWRD paper.  The
current implementation lives under `matlab/` and adds finite windows, matrix
estimators, derivative/primitive reconstruction, three-wire model families,
diagnostics, and measurement-degradation campaigns.

## Data Policy

No private measurement data are included in this repository.  The public tests
and campaigns generate synthetic signals locally.  Generated `.mat`, `.csv`,
figures, logs, and result folders are ignored by Git.

See [docs/DATA_POLICY.md](docs/DATA_POLICY.md) and
[docs/RELEASE_STRUCTURE.md](docs/RELEASE_STRUCTURE.md).

## Quick Start

From MATLAB:

```matlab
cd /path/to/SPACOR
startup_spacor
run_public_smoke
```

The smoke suite writes regenerated synthetic outputs under `results/`.  These
files are intentionally not tracked.

## Public API

```matlab
startup_spacor

% Synthetic single-phase fixture
[data, config, truth] = spacor.fixtures.singlephase_series_rl();
report = spacor.singlephase.characterize(data, config, struct('candidateMode','lti'));

% Synthetic three-wire fixture
[data3, config3, truth3] = spacor.fixtures.threewire_delta_g();
report3 = spacor.threewire.characterize(data3, config3, struct('candidateMode','mixed_blackbox'));
```

The main campaign entry points are:

```matlab
spacor.campaigns.run_singlephase_quick()
spacor.campaigns.run_threewire_quick()
spacor.campaigns.run_singlephase_degradation()
spacor.campaigns.run_threewire_degradation()
```

See [docs/API.md](docs/API.md) and [docs/CAMPAIGNS.md](docs/CAMPAIGNS.md).

## Citation

The original single-phase theoretical formulation is described in:

```bibtex
@article{montoya2021spacor,
  author  = {Montoya, Francisco G. and De Leon, Francisco and Arrabal-Campos, Francisco M. and Alcayde, Alfredo},
  title   = {Determination of Instantaneous Powers from a Novel Time-Domain Parameter Identification Method of Non-Linear Single-Phase Circuits},
  journal = {IEEE Transactions on Power Delivery},
  year    = {2021},
  doi     = {10.1109/TPWRD.2021.3133069}
}
```

Additional publications and the three-phase manuscript will be added to the
citation file as they are finalized.

## License

This fork preserves the upstream GPL-3.0 license.
