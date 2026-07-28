# NFC crystallography

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21639283.svg)](https://doi.org/10.5281/zenodo.21639283)

This repository is the ordinary, public entry point for a crystallographic
methods result developed in the NFC research program:

> Raw-only repeat-certified reciprocal-spot consolidation, multiscale
> primitive-lattice candidate generation, finite-index alternative testing,
> and explicit abstention.

The narrow empirical claim is that the method has real lattice-recovery
capability, including one correct prospective transfer on 8VTD, while also
having documented sensitivity limits and conservative misses. It is **not** a
claim that the complete D4–D7 pipeline transfers, an independent validation of
NFC, or evidence that NFC is a Theory of Everything.

## Scientific status

| Component | Current evidence |
| --- | --- |
| D4 spot signal | Supported on public real data; agrees strongly with conventional spot finding on 9JQ9. |
| D4.5 reciprocal consolidation | Causally useful on current controls; scan-local and increment-sign-aware. |
| Baseline B lattice method | Default comparator; real-data recovery capable with explicit abstention. |
| Candidate C phase substitution | Compatible continuous alternative; no prospective incremental value established. |
| Prospective record | 8VTD correct recovery; 9JZO conservative ambiguous decision with the physical family present. |
| D6c1–D7 transfer | Not established. |

The complete evidence table, including failures and exclusions, is generated
from `results/evidence_table.json`.

## Quick start

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[dev]'
pytest
nfc-cryst evidence --format markdown
```

Verify the method-source bindings:

```bash
nfc-cryst verify-methods
```

For a wheel installation, locate the bundled frozen release payload with:

```bash
nfc-cryst release-root
```

Run the repaired conventional qualification comparison on three DIALS
experiment files:

```bash
nfc-cryst qualify-conventional \
  --full full/refined.expt \
  --half-a half_a/refined.expt \
  --half-b half_b/refined.expt \
  --assigned-full 0.90 \
  --assigned-half-a 0.89 \
  --assigned-half-b 0.91 \
  --output qualification.json
```

The comparison minimizes orientation disagreement over a declared finite set
of proper unimodular reciprocal-basis transformations. It fails closed if no
metric-compatible transformation is found. This fixes the direct-matrix
comparison defect that excluded 6CKT and 6TPI; it does not alter baseline B.

## Reproducing a method decision

`methods/baseline_B/runtime/` preserves the exact evaluator source layout
bound by the archived baseline-B release. A generic wrapper runs that
evaluator on a truth-free compact D4.5 case:

```bash
python scripts/run_baseline_case.py \
  --case path/to/FIXED_D45_CASE.json.gz \
  --output result.json
```

Generating the case from public raw PILATUS MiniCBF frames requires a DIALS
runtime because the raw adapter consumes dxtbx geometry. The exact adapter and
frozen D4/D4.5 sources are retained under
`methods/baseline_B/raw_pipeline/` in a checkout and under the path printed by
`nfc-cryst release-root` in an installed wheel. Public data are not vendored.

The reproduction boundary is deliberate:

- compact evaluator and gate tests run in ordinary Python;
- raw-frame reconstruction requires DIALS/dxtbx;
- large public archives are downloaded separately and verified by manifest;
- historical releases remain the deep audit record.

See [REPRODUCING.md](REPRODUCING.md) for the full path.

## Release history

- `v0.1.0` is the immutable first public methods milestone. Its compact 8VTD
  replay remains valid, but the distributed package omitted frozen support
  modules needed for ordinary new-raw D4.5 construction.
- `v0.1.1` restores those exact hash-bound modules and verifies the installed
  wheel. It is a packaging-only correction and does not change baseline B,
  candidate C, any scientific threshold, or any historical outcome.

The version-independent archival identifier is
[Zenodo concept DOI 10.5281/zenodo.21639283](https://doi.org/10.5281/zenodo.21639283).

## Methods and comparators

Baseline B is controlling:

```text
ba18310a04a45c13f1fdf100599c77f2da9fa8ba2f43ec9e942719871b6edf48
```

Candidate C is an experimental, compatible alternative:

```text
e9e2b2872c64c37bc487824385488fb5b009c2739cddf5a4089409bebc39499d
```

An exact tie does not promote C. The paired cycle produced one valid tie on
9JZO and was closed after two further corpora exposed a representation defect
in the pre-execution conventional gate.

## Repository map

```text
src/nfc_cryst/               public utilities and invariant gate
methods/                     exact method bindings and frozen source
data_manifests/              public-input identities and roles
results/                     compact evidence and gate-development records
scripts/                     one-purpose reproduction commands
tests/                       fast scientific and integrity tests
paper/                       methods-paper draft
```

## Claim boundaries

- Development results remain development results even when they recover the
  deposited lattice.
- 8VTD is the current correct prospective recovery.
- 9JZO is a prospective conservative miss, not a wrong-lattice result.
- 6CKT and 6TPI were pre-execution gate exclusions, not method failures.
- Candidate C has no demonstrated incremental benefit.
- There has been no independent external reproduction yet.
- No crystallography result in this repository establishes NFC as a
  cosmology or Theory of Everything.

## Human–AI provenance

The research program and evidence packages were developed through an extended
human–AI collaboration led by Evan Thomas Kotler. Scientific claims here are
made on the basis of public inputs, executable code, preserved hashes, and
falsifiable outcomes—not on the authority of the conversations that produced
them. See [PROVENANCE.md](PROVENANCE.md).
