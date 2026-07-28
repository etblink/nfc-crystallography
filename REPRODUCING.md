# Reproducing the results

## Levels of reproduction

This repository separates three increasingly expensive checks.

### 1. Fast repository checks

```bash
python -m pip install -e '.[dev]'
pytest
nfc-cryst verify-methods
nfc-cryst evidence --format markdown
```

These checks validate the public utilities, method-source bindings, evidence
schema, and invariant conventional gate.

### 2. Compact truth-free D5 replay

Obtain a compact fixed-D4.5 case from one of the archived reproducibility
releases, then run:

```bash
python scripts/run_baseline_case.py \
  --case FIXED_D45_CASE.json.gz \
  --output baseline_result.json
```

The wrapper refuses a case whose truth basis is populated. It executes the
unchanged archived baseline evaluator and writes a canonical JSON result.

### 3. Public raw-frame replay

Raw-frame replay requires:

1. a public archive listed in `data_manifests/`;
2. a DIALS/dxtbx environment capable of reading the detector format;
3. the exact raw adapter in `methods/baseline_B/raw_pipeline/`;
4. the frozen D4 C decoder built locally;
5. sufficient disk, memory, and CPU time.

Public archives are deliberately not committed to Git. Use
`scripts/download_public_inputs.py` to download and verify a manifest that
contains a direct URL and SHA-256:

```bash
python scripts/download_public_inputs.py \
  data_manifests/8vtd.json \
  --destination downloads/8vtd.tar.bz2
```

The 8VTD and 9JZO manifests contain the whole public-archive SHA-256 values
observed during their completed acquisitions. A manifest without a whole-file
SHA-256 is not accepted by the convenience downloader; dataset-specific
source-file hashes remain preserved in the historical releases as an
additional layer.

## DIALS note

DIALS is not bundled or silently installed. Use an official DIALS release and
record its version. The repository gate reads DIALS `.expt` JSON directly and
does not require DIALS merely to compare already-refined crystal models.

## Historical audit artifacts

The exact archived ZIPs remain the controlling deep evidence:

- baseline B: `ba18310a…b6edf48`
- 8VTD prospective result: `de10c22c…c47b8`
- candidate C: `e9e2b287…c39499d`
- 9JZO paired result: `c5d43784…770202e`
- complete paired cycle: `e8135d3f…cf694d`

The repository does not rewrite those artifacts. It exposes their scientific
content in a form that can be run and reviewed without following the entire
historical governance chain.
