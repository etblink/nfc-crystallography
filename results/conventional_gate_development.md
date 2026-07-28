# Conventional qualification-gate development

The paired B/C cycle used direct orientation-matrix angles. That historical
gate excluded 6CKT and 6TPI before either NFC method ran:

| Corpus | Minimum assigned fraction | Maximum metric difference | Direct orientation difference |
| --- | ---: | ---: | ---: |
| 6CKT | 0.887896 | 0.006009 | 179.997951° |
| 6TPI | 0.916506 | 0.002438 | 179.998068° |

After the paired cycle closed, both corpora became gate-development inputs.
The replacement utility enumerates 3,480 proper unimodular reciprocal-basis
operators with entries in `[-1, 1]`, rejects operators outside the original
0.05 reciprocal-metric tolerance, and minimizes the physical polar-rotation
residual over the survivors.

| Corpus | Invariant maximum metric difference | Invariant maximum orientation difference | Development result |
| --- | ---: | ---: | --- |
| 6CKT | 0.005891 | 0.036427° | Stable |
| 6TPI | 0.002134 | 0.036397° | Stable |

The exact records are:

- `results/6ckt_invariant_gate_development.json`
- `results/6tpi_invariant_gate_development.json`

This does not retroactively convert either exclusion into a paired method
result. It establishes that the original exclusion mechanism was
representation-dependent and supplies a tested gate for a future prospective
cycle. The bounded operator family is part of the method record; invariance
beyond that declared family is not claimed.
