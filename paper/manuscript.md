# Raw-only reciprocal-spot consolidation and primitive-lattice recovery with explicit abstention

**Evan Thomas Kotler**

## Abstract

We report an exploratory crystallographic method that constructs
three-dimensional reciprocal observations from raw rotation images without
using a conventional unit cell or orientation, generates multiscale primitive
lattice candidates, audits finite-index alternatives, and can abstain when the
evidence does not identify a unique lattice. The pipeline has three stages:
D4 detects frame-local diffraction signal; D4.5 associates repeat detections
within acquisition series and consolidates them into reciprocal centroids; and
D5 generates, compares, and either accepts or rejects lattice families using
split reproducibility, held-out prediction, origin specificity, and
finite-index alternatives. D4 agreed strongly with conventional spot finding
on a public adverse corpus. D4.5 centroid consolidation causally repaired
real-data failures on two 6GN3 sweeps. The resulting D5 lineage recovered
physical primitive lattices on multiple development corpora and correctly
recovered 8VTD in a prospective run. A second prospective corpus, 9JZO,
produced an ambiguous decision despite containing the physical primitive
family among its surviving candidates. Other public corpora produced
abstentions or exposed an unrelated conventional-qualification defect. These
results establish real geometric recovery capability and explicit limitations;
they do not establish transfer of later intensity stages, independent external
validation, or support for NFC as a cosmological theory.

## 1. Introduction

Indexing diffraction images normally combines mature detector models, spot
finding, reciprocal-space mapping, lattice search, and refinement. This work
began from a different representation and asked a narrower empirical question:
can structure derived from raw frame-local patterns yield reproducible
primitive-lattice information on public macromolecular rotation datasets?

Early results on 9Z6F motivated a portability implementation for PILATUS
MiniCBF data. Open evaluation on 9JQ9 then produced an informative split
result. D4 detected genuine diffraction spots, including 58,528 mutual
one-pixel pairs with DIALS and a matched-intensity Spearman correlation of
0.903. The frozen D5 construction, however, was unstable, poorly separated
from rotation-preserving nulls, and forced lattice-like outputs on adverse
synthetic controls. Subsequent public controls separated two defects: repeated
frame-local detections required consolidation before lattice search, and the
lattice decision needed explicit alternative comparison and abstention.

The contribution reported here is therefore not a replacement for
conventional crystallographic processing. It is a raw-derived candidate and
evidence framework whose successful and unsuccessful cases can be tested
independently.

## 2. Methods

### 2.1 D4 frame-local signal

D4 consumes native detector pixels under a literal-pixel rule. The published
PILATUS adapter performs no resampling, cropping, or pixel-value scaling. It
maps detector and scan metadata through the native geometry and emits
frame-local peak observations. The exact CBF byte-offset decoder, D4 rule, and
adapter sources are hash-bound in the repository.

### 2.2 D4.5 repeat-certified consolidation

D4.5 operates independently within acquisition series. Traversal follows the
sign of the scan increment, associations do not cross scan-series boundaries,
and repeated detections are consolidated into signal-weighted reciprocal
centroids with deterministic ordering. The rule uses no conventional cell,
orientation, symmetry, or processed reflection object.

Causal ablation on 6GN3 distinguished consolidation from mere repeat
filtering. Unaggregated detections and certified original members failed,
whereas arithmetic and signal-weighted centroids recovered both full sweeps
and their independently constructed halves. A later scan-direction test found
that an ascending-only traversal produced no repeat objects for negative
increments; scan-local increment-sign-aware traversal repaired 9Z6F while
leaving positive-increment corpora unchanged.

### 2.3 Baseline B candidate and decision stack

Baseline B uses the archived multiscale and bidirectional successor with seed
floors 0.012, 0.009, and 0.006 Å⁻¹ and requires persistence at two scales.
The candidate generator and integer refinement are inherited from the frozen
D5 implementation. The successor adds:

1. independently fitted full, half-A, and half-B candidates;
2. primitive-lattice consensus across those fits;
3. held-out cross-half scoring;
4. zero-origin specificity against equal-density translated origins;
5. bidirectional index-two completion alternatives;
6. a complexity audit for finite-index aliases; and
7. `LATTICE_RECOVERED`, `AMBIGUOUS_LATTICE`, or
   `INSUFFICIENT_SIGNAL` outcomes.

Candidate C changes only the phase component to a fixed-scale Cauchy affinity
with scale 0.003 Å⁻¹. It preserved all current development classifications,
but a paired prospective run did not establish incremental value. Baseline B
therefore remains controlling.

### 2.4 Evidence roles

Corpora are classified by scientific exposure rather than a single binary
label. Development corpora may support mechanism discovery and ablation but
not out-of-sample transfer claims. A prospective run binds raw inputs, fixed
D4.5 feeds, method identity, and a truth-free decision before the conventional
cell and orientation are consulted for scoring. Pre-execution technical
exclusions are not counted as method outcomes.

### 2.5 Conventional qualification

The first paired cycle compared raw orientation matrices directly. This is not
invariant to equivalent reciprocal bases or indexing ambiguities and excluded
6CKT and 6TPI at nominal differences near 180° despite stable metrics and high
assigned fractions.

The repository implements a replacement qualification utility that searches a
declared finite set of proper unimodular integer basis transformations. For
each metric-admissible transform, it obtains the closest proper physical
rotation through polar decomposition and uses the minimum invariant
orientation residual. The utility is fail-closed and isolated from D4.5/D5.
It is a qualification-interface repair, not a change to baseline B.

## 3. Results

### 3.1 Signal and representation

D4’s agreement with conventional spot finding on 9JQ9 supports a genuine
detector-level signal. D4.5 centroid consolidation was the causal repair on
6GN3; repeat filtering alone was insufficient. Scan-local sign-aware traversal
recovered the historical 9Z6F primitive lattice across independent scan
partitions.

### 3.2 Recovery and abstention controls

The predecessor D5 could recover idealized positives but failed all six
structured abstention controls and could assign approximately 96% formal
support to a pure nonlattice cloud. The split-consensus and complexity-gated
successor recovered all then-established positive conditions and abstained on
adverse conditions in the development set. A later 6MFU miss exposed a fixed
seed-floor and one-sided finite-index limitation. Multiscale access plus a
bidirectional index-two audit resolved the alias on development controls
without introducing a confident wrong recovery.

These results establish capability, not universal validity. Several rules were
informed by the development corpora, and all such cases remain labeled
development evidence.

### 3.3 Prospective outcomes

Baseline B correctly recovered the physical primitive lattice on 8VTD after
its method and input identities were fixed and before conventional truth was
consulted. This is the strongest current positive transfer result.

On 9JZO, baseline B and candidate C returned the same
`AMBIGUOUS_LATTICE` decision with identical selected survivor families.
Post-decision scoring showed that the common latent family was
primitive-equivalent to the stable conventional solution. This was a
conservative missed recovery rather than a confident incorrect lattice.
Because the comparators tied, candidate C was not promoted.

6CKT and 6TPI never reached method execution. Both had high conventional
assignment and close reciprocal metrics, but failed a direct matrix-angle gate
at approximately 180°. They remain historical pre-execution exclusions and
are now used only to develop and test the invariant qualification interface.

### 3.4 Negative and unresolved cases

9JQ9 retained strong D4 signal but did not support a stable D5 lattice and
remains an unresolved adverse diagnostic. 5V0G and 6W61 produced
`INSUFFICIENT_SIGNAL` outcomes. Lower seed floors exposed near-reference metric
signal in these development cases but did not establish split-reproducible,
held-out lattice identification. Hard-cutoff and object-specific uncertainty
proxies were tested and rejected when they failed calibration or adverse
selectivity.

## 4. Discussion

The cumulative evidence supports a useful decomposition. D4 detects real
frame-local signal. D4.5 converts repeated detections into a more stable
three-dimensional reciprocal representation. The D5 candidate machinery can
reach physical primitive lattices on real public data. Final selection remains
the limiting inferential layer: a correct family may survive without earning a
recovery decision, while phase specificity alone is unsafe on geometry-error
and competing-lattice controls.

Explicit abstention is therefore a substantive scientific feature. A miss on a
stable positive control is a sensitivity limitation, but it is preferable to a
confident wrong lattice. The current prospective record contains one correct
recovery, one conservative miss with the physical family present, and no
confident prospective wrong recovery.

The results do not yet authorize conclusions about observation-level index
transfer, intensity prediction, or later D6c1–D7 stages. Nor do they provide a
quantitative basis for claims about NFC cosmology. The crystallographic method
must stand or fail as a reproducible algorithmic contribution.

## 5. Availability and next tests

The repository includes exact baseline source bindings, a compact truth-free
runner, public data manifests, fast tests, a complete evidence table, and the
invariant conventional gate. Large public raw archives are downloaded
separately. Historical ZIP releases remain the deep audit record.

The next empirical test should use one untouched public corpus, qualify it
with the invariant conventional gate, and execute unchanged baseline B before
truth scoring. Exploratory D6c1 work should begin only after another correct
prospective geometric recovery on a corpus that also supports a reproducible
public observation-level reduction.
