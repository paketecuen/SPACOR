# Public MATLAB API

This document describes the public, synthetic-data API exposed by the current
repository.  Private real-data loaders are not included in this release.

## Startup

```matlab
cd /path/to/SPACOR
startup_spacor
```

This adds `matlab/`, `scripts/`, and `tests/` to the MATLAB path.

## Single-Phase Identification

```matlab
[data, config, truth] = spacor.fixtures.singlephase_series_rl();
report = spacor.singlephase.characterize(data, config, ...
    struct('candidateMode','lti'));
```

The returned `report` contains the selected equivalent, candidate ranking,
window diagnostics, residuals, and recovered physical parameters.

## Three-Wire Identification

```matlab
[data, config, truth] = spacor.fixtures.threewire_delta_g();
report = spacor.threewire.characterize(data, config, ...
    struct('candidateMode','mixed_blackbox'));
```

The three-wire API evaluates candidate delta and wye equivalent families from
terminal line-to-line voltages and line currents.  The implementation reports
rank, condition number, KCL/KVL consistency, terminal-power residuals, energy
residuals, passivity checks, and physical parameter estimates.

## Fixtures

Synthetic fixtures are generated locally:

```matlab
manifest = spacor.fixtures.write_synthetic_fixtures();
```

The generated `.mat` files are useful for local tests, but they are not tracked
in Git.

## Campaigns

```matlab
spacor.campaigns.run_singlephase_quick()
spacor.campaigns.run_threewire_quick()
spacor.campaigns.run_singlephase_degradation()
spacor.campaigns.run_threewire_degradation()
```

Campaign outputs are written under `results/` and ignored by Git.

## Tests

```matlab
run_public_smoke
```

The public smoke suite uses synthetic fixtures and synthetic degradation
campaigns only.
