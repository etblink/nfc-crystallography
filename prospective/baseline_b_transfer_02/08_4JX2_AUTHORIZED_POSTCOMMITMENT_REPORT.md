# 4JX2 authorized post-commitment scoring report

## Immutable bindings verified

| Binding | Verified value |
| --- | --- |
| Commitment file SHA-256 | `7e9e93a0a3c8c77708b6820f8dc064c2fe66fa9d7dedf808cb74f14023a1ed1a` |
| Commitment semantic SHA-256 | `2aac64f387f036f887146dc0fb1b5603dba680c1a18fbedf76b87be720fd6609` |
| Bound truth-free evaluation SHA-256 | `83f927a9311504184484eacb9ae1671876f38c16681cd106e3d7954b5d9d129c` |
| Committed baseline-B decision | `INSUFFICIENT_SIGNAL` |
| Committed reason | `NO_PERSISTENT_FAMILY` |
| Candidate C execution | `NOT RUN` |

The commitment and the already-recorded truth-free result were read and hash-verified before conventional truth was used. Baseline B was not rerun, Candidate C was not run, and the committed decision remains unchanged.

## Conventional full/split qualification evidence

All three independent cell-blind DIALS indexings produced one crystal model. The solution values below are from the sealed conventional artifacts.

| Partition | DIALS experiments | Hall symbol | Unit cell (Å, °) | Assigned fraction |
| --- | ---: | --- | --- | ---: |
| FULL | 11 | `P 1` | 63.7600, 121.5253, 144.6206, 90.0836, 90.0765, 90.0046 | 0.735408 |
| HALF_A | 6 | `P 1` | 63.7488, 121.4742, 144.5888, 90.1006, 89.9201, 89.9943 | 0.733904 |
| HALF_B | 6 | `P 1` | 63.7733, 121.5565, 144.6409, 90.0913, 89.9240, 89.9826 | 0.737539 |

The released basis-invariant conventional gate passed. Its fixed limits and observed worst-case values were:

| Gate quantity | Limit | Observed |
| --- | ---: | ---: |
| Assigned fraction, each partition | ≥ 0.500000 | ≥ 0.733904 |
| Relative metric difference | ≤ 0.050000 | 0.00083646 |
| Basis-invariant orientation difference | ≤ 2.000000° | 0.12206669° |
| Finite basis-search bound | 1 | 1 |

The frozen gate outcome was `CONVENTIONAL_FULL_AND_SPLIT_INDEXABILITY_STABLE`.

The reciprocal orientation matrices (Å⁻¹, DIALS column-basis convention) were:

```text
FULL
[-0.013134, -0.002109,  0.003328]
[-0.006371, -0.002463, -0.005977]
[ 0.005736, -0.007563,  0.001007]

HALF_A
[-0.013140,  0.002109, -0.003325]
[-0.006370,  0.002459,  0.005981]
[ 0.005729,  0.007568, -0.001003]

HALF_B
[-0.013124,  0.002107, -0.003334]
[-0.006379,  0.002471,  0.005972]
[ 0.005739,  0.007558, -0.001011]
```

## Primitive-equivalence diagnostic only

The frozen baseline primitive-equivalence routine was used without modification: maximum transform-integer deviation `0.03`, required absolute determinant `1`. Each comparison is matched to the same FULL/HALF_A/HALF_B partition. This is post-commitment diagnostic scoring; it is not a baseline-B decision input and does not promote any latent family.

| Already-recorded item | Scale (Å⁻¹) | FULL (deviation, det) | HALF_A (deviation, det) | HALF_B (deviation, det) | Primitive-equivalent across all three? |
| --- | ---: | --- | --- | --- | --- |
| Direct candidate | 0.012 | 0.481530, 1 | 0.484421, 1 | 0.481938, 1 | No |
| Direct candidate | 0.009 | 0.488455, 1 | 0.494271, 2 | 0.492258, 2 | No |
| Direct candidate | 0.006 | 0.014181, 1 | 0.014835, 1 | 0.012388, 1 | Yes |
| Weak completion family | 0.006 | 0.014107, 2 | 0.014764, 2 | 0.012330, 2 | No |

No other admitted, surviving, weak-conflict, strong-single-scale-conflict, or direct-single-scale-conflict family was recorded. The matching 0.006 direct candidate is a **single-scale latent family only**. It did not meet the frozen persistence rule and therefore remains non-promoted.

## Final classification

`CONSERVATIVE_MISSED_RECOVERY`

The prospective baseline abstained (`INSUFFICIENT_SIGNAL`) rather than making an incorrect recovery. The post-commitment diagnostic shows that its already-recorded single-scale 0.006 Å⁻¹ direct candidate was primitive-equivalent to the qualifying conventional FULL, HALF_A, and HALF_B solutions. This supports the missed-recovery classification, but does **not** change the baseline decision to `LATTICE_RECOVERED`.

## Observation-level reduction and D6c1

No public reproducible observation-level reduction exists. The locally sealed conventional-processing products and this post-commitment comparison do not constitute that separate public reduction. Accordingly, the D6c1 gate does not change and remains unsatisfied.

## Reproducibility record

The machine-readable read-only diagnostic record is `4JX2_authorized_postcommit_primitive_equivalence_diagnostic.json`.

| Artifact | SHA-256 |
| --- | --- |
| Diagnostic JSON | `2df40d15d947c0a5898f67594a214d1a98948bc25a0674e72ee8c2f1e56442b8` |
| Frozen equivalence wrapper | `2fe1c017d287ee7cb24abb8c8c7cc100f46c72c24f5b38b6ad0654b456a06002` |
| Frozen equivalence delegate | `0f142754f11e0703642a304f74452dde9233ffb4073709afe206b8d07ddec2f5` |
