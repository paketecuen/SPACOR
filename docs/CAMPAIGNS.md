# Public Campaigns

The public campaign layer is designed to regenerate synthetic evidence without
shipping measurement data.

## Quick Campaigns

```matlab
startup_spacor
spacor.campaigns.run_singlephase_quick()
spacor.campaigns.run_threewire_quick()
```

Quick campaigns verify the nominal behavior of the single-phase and three-wire
identification APIs on deterministic synthetic waveforms.

## Measurement-Degradation Campaigns

```matlab
startup_spacor
spacor.campaigns.run_singlephase_degradation()
spacor.campaigns.run_threewire_degradation()
```

The degradation campaigns perturb synthetic terminal waveforms with controlled
measurement effects such as noise, quantization, gain/offset errors, window
length variation, and filtering/delay effects where implemented by the runtime.

## Output Files

Campaigns write regenerated outputs under `results/`, typically including:

```text
summary.csv
parameters.csv
stability.csv
information.csv
method_ranking.csv
candidates.csv
report.mat
figures/
```

These outputs are reproducible artifacts and are not tracked by Git.

## Private Measurements

Private laboratory or field measurements are intentionally absent from this
release.  Future public data packages should be added as explicit, documented
datasets with consent, provenance, and separate checksums.
