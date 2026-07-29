# Raw-only repeat-certified reciprocal-spot consolidation and multiscale primitive-lattice recovery with explicit abstention

**Evan Thomas Kotler**

Independent researcher, Las Vegas, Nevada, USA

## Abstract

Indexing macromolecular rotation data ordinarily combines detector modeling,
spot finding, reciprocal-space mapping, lattice search, and refinement. We
evaluate an exploratory alternative that derives three-dimensional reciprocal
observations from raw images without using a conventional unit cell,
orientation matrix, symmetry assignment, or processed reflection file. The
pipeline has three stages. D4 detects frame-local diffraction signal using a
frozen robust pixel rule. D4.5 associates repeat detections only across
consecutive frames within an acquisition series and consolidates them into
signal-weighted reciprocal centroids. D5 generates primitive-lattice candidates
at three reciprocal-length floors and either recovers a persistent family or
abstains after split-reproducibility, held-out origin-specificity, finite-index,
and complexity tests.

On 9JQ9, D4 produced 58,528 mutual one-pixel matches with conventional DIALS
spots, with a matched-intensity Spearman correlation of 0.903. Controlled
ablations on two 6GN3 sweeps showed that centroid consolidation, rather than
repeat filtering alone, caused recovery by the unchanged downstream lattice
method. The developed D5 lineage recovered primitive lattices on several public
development corpora and passed its declared synthetic recovery and adverse
controls without a confident incorrect recovery. In prospective use, the
released baseline correctly recovered 8VTD. It abstained on 9JZO and 4JX2,
although post-commitment scoring found the physical primitive family among
9JZO's latent survivors and at one of three tested scales for 4JX2. These remain
missed recoveries, not retroactive successes. Across three prospective
corpora, the record is therefore one correct recovery, two conservative misses,
and no confident wrong recovery.

The results support a raw-derived reciprocal representation and a
recovery-or-abstention framework with real transfer capability and substantial
sensitivity limits. They do not establish broad indexing superiority,
observation-level or intensity transfer, independent validation, or any
cosmological claim.

**Keywords:** diffraction indexing; reciprocal lattice; rotation data;
spot finding; abstention; reproducible research

## 1. Introduction

Indexing a rotation diffraction experiment requires a consistent explanation
of observed spot positions by a reciprocal lattice and an orientation. Mature
programs such as XDS and DIALS combine detector geometry, spot finding, lattice
search, indexing, and refinement in robust processing systems [1,2]. Established
indexing approaches include Fourier analysis of reciprocal-space difference
vectors [3], robust real-space grid searches [4], and extensions for multiple
lattices [5]. Reduced-cell and basis-equivalence machinery is essential when
comparing solutions represented by different valid lattice bases [6].

The present work began with a narrower question: can a raw-image
representation developed independently of a conventional indexing solution
produce reproducible primitive-lattice information on public macromolecular
rotation datasets? The question is empirical rather than terminological. A
useful result must recover physical lattices outside its development examples,
reject or abstain on unsupported inputs, expose its failures, and be
reproducible without accepting the broader theoretical vocabulary that
motivated the investigation.

Initial work on 9Z6F used a detector-specific implementation. A portability
layer was subsequently added for native PILATUS 1.2 MiniCBF frames, retaining
native detector dimensions and geometry and prohibiting resampling, cropping,
or pixel-value scaling. PILATUS detectors and the CBF/imgCIF family are
well-established crystallographic technologies [7,8]. Public raw datasets were
obtained from the Integrated Resource for Reproducibility in Macromolecular
Crystallography (IRRMC) and associated ProteinDiffraction repositories [9].

Open evaluation on 9JQ9 separated a positive upstream result from a negative
downstream result. Frame-local D4 detections agreed strongly with conventional
spots, but the original D5 lattice objective was unstable and could assign high
formal support to nonlattice controls. Further public and synthetic controls
isolated three problems: repeated frame-local detections required
three-dimensional consolidation; a fixed reciprocal-length floor could exclude
large-cell modes; and a forced-recovery objective required explicit alternative
comparison and abstention.

This paper reports the resulting public method, called **baseline B**, and the
complete evidence balance supporting it. Baseline B comprises a fixed D4.5
representation, multiscale candidate access, bidirectional index-two
alternatives, split and held-out tests, a complexity audit, and three possible
decisions: `LATTICE_RECOVERED`, `AMBIGUOUS_LATTICE`, or
`INSUFFICIENT_SIGNAL`. The contribution is not presented as a production
replacement for conventional indexing. It is a raw-derived candidate and
evidence framework with one correct prospective transfer, two prospective
misses, and documented limits.

## 2. Materials and methods

### 2.1 Study design and evidence roles

All scientific roles were assigned according to information exposure. A
development corpus could be used to diagnose or change the method and therefore
could not later support an out-of-sample claim for that version. A prospective
corpus was selected before method execution; raw inputs, fixed method identity,
truth-free inputs, and the decision were bound before the conventional unit
cell and orientation were consulted for scoring. A pre-execution exclusion did
not count as a method outcome.

The evidence roles are summarized in Table 1. Separate sweeps of 6GN3 were
treated as distinct experimental conditions but not as distinct underlying
structures. The small, deliberately heterogeneous corpus collection is an
evaluation series, not a population sample from which a general indexing
success rate can be estimated.

**Table 1. Public real-data evidence roles and final interpretations.**

| Corpus or condition | Role at the reported stage | Committed or controlling outcome | Interpretation |
| --- | --- | --- | --- |
| 9Z6F | Method development | `LATTICE_RECOVERED` | Historical development recovery |
| 9JQ9 | Adverse exploratory diagnostic | D4 signal; D5 not established | Real spot signal; unstable lattice evidence and unresolved provenance |
| 6GN2 | Development positive control | `LATTICE_RECOVERED` | Public real-data primitive-lattice recovery |
| 6GN3 sweep 1 | Development positive control and ablation | `LATTICE_RECOVERED` after D4.5 consolidation | Representation repair |
| 6GN3 sweep 2 | Development positive control and ablation | `LATTICE_RECOVERED` after D4.5 consolidation | Representation repair |
| 6MFU | Prospective then development | Abstention; later post-hoc recovery | Seed-floor and index-two alias diagnostic |
| 8VTD | Prospective out-of-sample | Correct `LATTICE_RECOVERED` | Correct prospective primitive recovery |
| 9JZO | Paired prospective | `AMBIGUOUS_LATTICE` for B and C | Conservative miss; latent primitive family present |
| 4JX2 | Baseline-B-only prospective | `INSUFFICIENT_SIGNAL` | Conservative miss; primitive family present at one tested scale |
| 5V0G | Truth-exposed development | `INSUFFICIENT_SIGNAL` | Known-positive development miss |
| 6W61 | Staged multisweep prospective then development | `INSUFFICIENT_SIGNAL` | Archived conservative abstention |
| 4G2A | Pre-execution exclusion | No NFC execution | Unsupported SLS 1.0 / X-CW raw interface |
| 6CKT, 6TPI | Pre-execution exclusions | No NFC execution | Historical qualification-gate defect; later gate-development corpora |

### 2.2 D4 frame-local detection

D4 consumes decoded native detector values. The public portability path
supports signed 32-bit `x-CBF_BYTE_OFFSET` PILATUS MiniCBF data through dxtbx
geometry and a byte-exact decoder. Masked or invalid pixels are excluded.
Native detector dimensions, pixel size, distance, beam center, axes,
wavelength, and scan angle are retained; no resampling, cropping, or
pixel-value scaling is performed.

Each image is divided into 32 × 32 pixel tiles. For the nonnegative values in
tile \(t\), the background \(m_t\) is the median and the robust scale is

\[
s_t=\max\left\{1,\;1.4826\,\mathrm{MAD}_t,\;
\sqrt{\max(m_t,0)+1}\right\}.
\]

For pixel value \(I_i\) in tile \(t(i)\),

\[
z_i=\frac{I_i-m_{t(i)}}{s_{t(i)}}.
\]

A connected component begins from a pixel satisfying \(z_i\geq 8\) and
\(I_i-m_{t(i)}\geq 8\), and grows through eight-connected pixels satisfying
\(z_i\geq 3\) and positive excess. A component must contain at least two
pixels and have integrated positive excess divided by formal noise of at least
12. Multiple-local-maximum components remain typed as unresolved rather than
being split post hoc. The peak centroid is the positive-excess-weighted pixel
center. The component footprint, half-pixel extent, and formal centroid terms
define a deterministic enclosure used for association; this enclosure is not
claimed to be a calibrated confidence region.

### 2.3 Reciprocal-coordinate mapping

Let \(\mathbf f\) and \(\mathbf s\) denote the native detector fast and slow
unit vectors, \(\mathbf n\) the detector normal, \(D\) the detector distance,
\((b_x,b_y)\) the beam center in pixel coordinates, \(p_x,p_y\) the native
pixel dimensions, and \((x,y)\) the zero-based D4 centroid. The detector point
used for ray construction is

\[
\mathbf d=D\mathbf n+
p_x(x+0.5-b_x)\mathbf f+
p_y(y+0.5-b_y)\mathbf s .
\]

For wavelength \(\lambda\), incident-beam direction \(\hat{\mathbf b}\), and
sample orientation \(\mathbf R_\phi\) at the frame midpoint,

\[
\mathbf s_0=\frac{\hat{\mathbf b}}{\lambda},\qquad
\mathbf s_1=\frac{\mathbf d/\|\mathbf d\|}{\lambda},\qquad
\mathbf q=\mathbf R_\phi^{\mathsf T}(\mathbf s_1-\mathbf s_0).
\]

Reciprocal coordinates are expressed in cycles Å\(^{-1}\), without a
\(2\pi\) factor. The mapping uses the scan-axis direction and sign reported by
the image model. A public SMV projection qualification on 501 frames compared
56,491 coordinates with a reference implementation, with no array or
source-hash mismatch and a maximum component difference of
\(2.29\times10^{-15}\) Å\(^{-1}\).

### 2.4 D4.5 repeat-certified consolidation

D4.5 is constructed independently for the full scan and for chronological
half-A and half-B feeds. Processing is local to an acquisition series, follows
the sign of the scan increment, and prohibits cross-series association.
Candidate associations are considered only between consecutive frames. Two
detections are compatible when their reciprocal enclosures overlap. Ambiguous
many-candidate nodes are rejected; remaining one-to-one links form
chronological repeat paths.

For repeat path \(P\), reciprocal centroid \(\bar{\mathbf q}_P\) is

\[
\bar{\mathbf q}_P=
\frac{\sum_{i\in P} w_i\mathbf q_i}{\sum_{i\in P}w_i},
\]

where \(w_i\) is the D4 integrated positive signal. An aggregate is admitted
only if its centroid remains inside every member's formal reciprocal
enclosure and the required contiguous partition is unique. Output ordering is
deterministic. The construction uses no conventional unit cell, orientation,
symmetry, processed reflection file, or truth-derived filter.

### 2.5 Multiscale D5 candidate generation

The D5 candidate kernel is inherited unchanged from the frozen exploratory
implementation. It balances at most 1,000 reciprocal objects per scan,
collects 15 nearest-neighbor difference vectors per object up to
0.12 Å\(^{-1}\), bins the differences at a width of 0.003 Å\(^{-1}\), retains
local modes at 10% of the maximum bin count, deduplicates modes within
0.004 Å\(^{-1}\), and forms independent triples subject to minimum pair-sine
0.1 and minimum normalized determinant 0.05. Integer refinement uses a fixed
descending residual schedule from 0.014 to 0.003 Å\(^{-1}\).

Baseline B repeats this search with lower difference-vector floors

\[
0.012,\quad 0.009,\quad 0.006\;\text{Å}^{-1},
\]

which equal four, three, and two times the existing voxel width. Each scale is
fitted independently on the full, half-A, and half-B feeds. A candidate family
must occur at at least two scales to authorize recovery. A single-scale family
may expose a conflict but can never produce `LATTICE_RECOVERED`.

### 2.6 Split, phase, finite-index, and complexity evidence

Primitive consensus is evaluated on the full and both half-scan bases. For
bases \(\mathbf B_1\) and \(\mathbf B_2\), the nearest integer transform to
\(\mathbf B_1^{-1}\mathbf B_2\) must have absolute determinant one and maximum
elementwise deviation no greater than 0.03 for all three pairwise
comparisons.

For a basis \(\mathbf B\), reciprocal observation \(\mathbf q\), and
fractional origin \(\mathbf o\), the phase residual is

\[
r(\mathbf q;\mathbf B,\mathbf o)
=\min_{\mathbf h\in\mathbb Z^3}
\|\mathbf q-\mathbf B(\mathbf h+\mathbf o)\|_2 .
\]

Baseline B measures the fraction of observations with
\(r\leq0.003\) Å\(^{-1}\). The zero origin is compared with the 63 nonzero
origins in
\(\{0,\tfrac14,\tfrac12,\tfrac34\}^3\). Every null uses the identical basis
and reciprocal-grid density. The zero-origin support must exceed the null
median by at least 0.1 and have empirical rank \(p\leq0.03125\). These tests
are required for the full basis on the full feed, each half basis on its own
feed, and both cross-half held-out directions.

The bidirectional audit evaluates the seven index-two primitive-completion
branches \(\mathbf B\mathbf H^{-1}\), complementing coarsening checks of
\(\mathbf B\mathbf H\). Paired support must satisfy \(p\leq0.03125\).
Completion evidence is classified by fixed minimum gains of 0.1 (strong) or
0.015625 (weak), with the same weak branch required across all five own and
held-out comparisons. Separately, all three-dimensional upper Hermite normal
form coarsenings with integer index 2 through 8 are tested. A lower-complexity
countermodel that retains at least 95% of the selected support and passes the
same phase test prevents recovery.

Decision precedence is fixed: insufficient objects or candidate failure;
full/split lattice inconsistency; phase failure; supported lower-complexity or
finite-index ambiguity; otherwise recovery. The only allowed final labels are
`LATTICE_RECOVERED`, `AMBIGUOUS_LATTICE`, and `INSUFFICIENT_SIGNAL`.

### 2.7 Comparator C

Candidate C changes only the phase statistic to a continuous Cauchy affinity
with fixed scale 0.003 Å\(^{-1}\). It leaves D4.5, candidate generation,
floors, persistence, split rules, finite-index alternatives, complexity
checks, thresholds, and abstention logic unchanged. C reproduced all 43
development-set classifications and all 17 selected recovery families of B.
A prospective paired test on 9JZO produced identical final and internal family
decisions. Because no incremental value was demonstrated, baseline B remains
the controlling method and C is reported only as an experimental compatible
alternative.

### 2.8 Conventional qualification and post-decision scoring

Conventional indexing was used to establish that a proposed prospective
corpus was an interpretable positive control, not to generate or filter NFC
candidates. Qualification required stable cell-blind full and split indexing,
assigned fractions of at least 0.5, pairwise reciprocal-metric disagreement
no greater than 0.05, and orientation disagreement no greater than 2°.

The current comparison is invariant over a declared finite family of proper
unimodular reciprocal-basis changes. For every metric-admissible transform it
uses polar decomposition to obtain the closest proper physical rotation and
reports the minimum angular residual. With basis-search entries restricted to
\(\{-1,0,1\}\), 3,480 proper unimodular transforms are examined. The
comparison fails closed if no admissible transform exists. This utility is
isolated from baseline B and was developed after the historical exclusions of
6CKT and 6TPI; those exclusions were not rewritten.

After a truth-free decision was cryptographically bound, candidate and
conventional bases were compared modulo an integer transformation. A
primitive-equivalent score required a nearest integer transform of absolute
determinant one within the frozen numerical tolerance. A latent or
single-scale primitive-equivalent family was recorded diagnostically but did
not change an abstention into a recovery.

### 2.9 Implementation and reproducibility

The public implementation is Python 3.10 or later with NumPy and SciPy; the
raw-frame path additionally requires a lawful DIALS/dxtbx runtime [2]. The
released baseline archive has SHA-256
`ba18310a04a45c13f1fdf100599c77f2da9fa8ba2f43ec9e942719871b6edf48`,
and every bundled source used
by the method is individually hash-bound. Public release `v0.1.1` is a
packaging-only correction to `v0.1.0`: it restores two exact frozen support
modules omitted from the first distributed package and changes no scientific
rule or result.

Fast repository tests verify source bindings, evidence records, the invariant
gate, and an installed-package raw-builder import. The compact 8VTD replay has
zero nonfloating mismatches against its archived decision and maximum absolute
floating disagreement \(2.83\times10^{-15}\). Large public raw archives are
not vendored; their URLs, byte sizes, and SHA-256 identities are recorded in
machine-readable manifests.

## 3. Results

### 3.1 D4 detects real frame-local diffraction signal

On the substantially different 9JQ9 detector and specimen, D4 produced 58,528
mutual pairs within one pixel of cell-blind DIALS spot centroids. Matched
integrated intensities had Spearman correlation 0.903. This establishes that
D4 responds to real diffraction spots, but it does not show superiority to
DIALS or incremental scientific information beyond conventional spot finding.

### 3.2 Centroid consolidation is a causal representation repair

The 6GN3 ablation held the downstream frozen D5 method fixed while varying the
input representation. For both sweeps, unaggregated D4 detections failed.
Restricting the feed to members of certified repeat paths without moving the
points also failed. Arithmetic and signal-weighted repeat centroids recovered
the physical primitive lattice on both full sweeps and their independently
constructed halves. Thus repeat filtering alone was insufficient, whereas
three-dimensional centroid consolidation was the supported causal repair.
Signal weighting was not necessary for the observed recovery, although it is
the fixed public rule.

An independent portability test isolated an ascending-scan assumption:
ascending-only traversal returned no repeat objects on negative-increment
9Z6F. Scan-local, increment-sign-aware traversal produced 17,880 repeat
objects and recovered the historical primitive lattice on the full corpus and
both scan-configuration partitions, while preserving the existing positive-
increment aggregate arrays.

### 3.3 Development controls motivate recovery-or-abstention

The predecessor objective could recover four idealized synthetic positive
cases but returned no required structured abstention in six adverse cases. In
a pure nonlattice cloud it retained 1,767 weak modes and assigned 95.625%
formal support, demonstrating that reciprocal-grid density could manufacture
a high support fraction.

The developed multiscale, bidirectional stack produced the required outcome
on all 32 cases in its fixed development census: 16
`LATTICE_RECOVERED`, two `AMBIGUOUS_LATTICE`, and 14
`INSUFFICIENT_SIGNAL`, with no incorrect recovery. The six targeted alias
controls included a balanced large cell, pseudotranslation, a true index-two
sublattice, a sparse superlattice coset, equal competing completions, and a
large-cell nonlattice background. Recovery was required for the first three;
structured abstention was required for the latter three.

These are development-set results, not an unbiased estimate of transfer.
Their scientific value is mechanistic: multiscale access repaired a
demonstrated floor exclusion, and the bidirectional audit prevented the
lowered floors from accepting a sparse parent explanation.

### 3.4 Real-data development recoveries and adverse cases

The method lineage recovered primitive lattices on 9Z6F, 6GN2, both 6GN3
sweeps, and 6MFU after 6MFU entered development. The 6MFU diagnosis was
specific: a required primitive-completion vector of 0.0111337 Å\(^{-1}\) lay
below the original 0.012 Å\(^{-1}\) floor, and the physical completion lay in
the direction \(\mathbf B\mathbf H^{-1}\) that a coarsening-only audit did not
search. The original prospective 6MFU abstention remains a miss.

9JQ9 remained adverse after D4.5 consolidation. Although it supplied 32,676
repeat-certified objects, downstream fits were wrong and split-unstable and
the public raw-to-deposited provenance could not be fully resolved. 5V0G and
6W61 also returned `INSUFFICIENT_SIGNAL`. Lower floors made a near-reference
full-scan metric visible in 5V0G and exposed physical-metric signal in parts of
6W61, but did not establish the split-reproducible, held-out evidence required
for recovery. These corpora were used to localize sensitivity, not to
reinterpret the archived decisions.

### 3.5 Prospective outcomes

Table 2 gives the complete prospective record for the public baseline after
development. The decision, not the presence of a plausible latent family,
determines whether a run counts as a recovery.

**Table 2. Prospective real-data outcomes.**

| Corpus | Method execution | Committed result | Post-commitment assessment |
| --- | --- | --- | --- |
| 8VTD | Baseline B | `LATTICE_RECOVERED` | Correct physical primitive lattice |
| 9JZO | Baseline B and candidate C | Both `AMBIGUOUS_LATTICE` | Exact comparator tie; common latent family primitive-equivalent |
| 4JX2 | Unchanged public `v0.1.1` baseline B | `INSUFFICIENT_SIGNAL` / `NO_PERSISTENT_FAMILY` | 0.006 Å\(^{-1}\) direct family primitive-equivalent; present at one scale only |

The 8VTD decision is the strongest positive transfer result. Its method and
input identities were fixed before conventional truth was consulted, and the
committed family was primitive-equivalent to the conventional solution.

On 9JZO, B and C produced identical `AMBIGUOUS_LATTICE` decisions and
identical survivor families. Post-commitment scoring found their common latent
family primitive-equivalent. The result is a conservative miss and supplies
no incremental evidence for C.

The 4JX2 archive contained 542 native CBF frames; 540 entered the fixed feeds.
D4.5 produced 31,668 full-feed objects, 16,024 half-A objects, and 15,647
half-B objects. Baseline B returned `INSUFFICIENT_SIGNAL` with reason
`NO_PERSISTENT_FAMILY`. Per-scale decisions were `INSUFFICIENT_SIGNAL` at
0.012 Å\(^{-1}\), `AMBIGUOUS_LATTICE` at 0.009 Å\(^{-1}\), and direct
recovery at 0.006 Å\(^{-1}\). The truth-free commitment was bound before
reveal (semantic SHA-256
`2aac64f387f036f887146dc0fb1b5603dba680c1a18fbedf76b87be720fd6609`).
Post-commitment scoring found the 0.006 Å\(^{-1}\) direct bases
primitive-equivalent to the conventional full and both split solutions, with
maximum transform deviations 0.01418, 0.01483, and 0.01239. The family was
not present at two scales, so the committed abstention is preserved.

The aggregate prospective record is one correct recovery, two conservative
misses, and zero confident wrong recoveries. With only three corpora, it does
not support a general success-rate estimate.

### 3.6 Conventional-gate repair

4G2A was excluded even earlier in the same intake sequence. All 578 raw
frames used an `SLS_1.0` header and `X, CW` rotation-axis convention not
admitted by the released fixed raw interface. No D4.5 feed or NFC decision was
produced, so this is an interface exclusion rather than a lattice-method
outcome.

6CKT and 6TPI were excluded before NFC execution because a direct orientation
matrix comparison reported differences near 180°. After the paired cycle was
closed, the basis-invariant utility found maximum full/split orientation
differences of 0.03643° and 0.03640°, respectively, with minimum assigned
fractions 0.8879 and 0.9165. This supports a representation-dependent defect
in the original qualification gate. The repair prevents future equivalent
indexing representations from failing for that reason; it creates no
retroactive NFC result for either corpus.

## 4. Discussion

### 4.1 Supported scientific claims

The cumulative evidence supports three distinct claims. First, D4 detects real
frame-local diffraction signal on public data. Second, D4.5 repeat-centroid
consolidation can convert fragmented frame-local detections into a reciprocal
representation useful for lattice search without importing a unit cell or
orientation. Third, the multiscale candidate machinery can reach physical
primitive lattices on real corpora, including one prospective transfer.

The results also support explicit abstention as an essential part of the
method. The predecessor's high support on a pure nonlattice control showed
that forced output and reciprocal-grid density can create false confidence.
Split reproducibility, held-out origin tests, bidirectional finite-index
alternatives, and complexity countermodels reject different failure modes and
are not redundant decorations around one score.

The prospective misses identify a sensitivity boundary rather than a
false-recovery boundary. In both 9JZO and 4JX2, candidate generation reached a
physical primitive family, but the final evidence rules did not authorize
recovery. That fact is useful for future successor research, but it cannot be
counted as present-method success. Conversely, no confident wrong prospective
lattice has yet been observed. Abstention on a stable positive control is a
scientific miss, even when safer than returning the wrong answer.

### 4.2 Limitations

The evaluation has important limitations.

1. The public raw path currently supports one new detector family, native
   PILATUS 1.2 MiniCBF, rather than the breadth of detector formats supported
   by established crystallographic suites.
2. The three-corpus prospective record is too small to estimate general
   recovery, error, or abstention rates.
3. Development recoveries are concentrated in a small number of underlying
   structures, and several thresholds and alternative families were informed
   by those controls.
4. D4 was compared closely with conventional spot finding on 9JQ9, but the
   present work is not a broad accuracy, speed, or resource benchmark against
   DIALS, XDS, or other indexers.
5. D4's formal reciprocal enclosure is deterministic, not calibrated
   object-specific uncertainty. Attempts to use enriched formal scale directly
   as a likelihood width did not improve held-out discrimination.
6. The origin-null family is finite; finite-index completions are restricted
   to index two; and the complexity audit tests upper-Hermite-normal-form
   coarsenings only through index eight. Candidate generation does not
   enumerate every possible lattice.
7. Baseline B contains hard floors and thresholds. It missed physical families
   that were latent or visible at only one tested scale.
8. Large raw archives must be downloaded from external public repositories;
   the software release is therefore reproducible from public sources but not
   fully self-contained.
9. No independent external clean-install reproduction has yet been reported.
10. No result here establishes observation-level index transfer, intensity
    prediction, anomalous handling, symmetry assignment, structure
    determination, D6c1-D7 transfer, or a broader NFC physical theory.

### 4.3 Implications and next tests

Baseline B should remain fixed as the reported method. Candidate C is a
mathematically continuous phase alternative but has no demonstrated
incremental value. Further tuning on 9JZO, 4JX2, 5V0G, or 6W61 would convert
those corpora into development inputs without supplying new evidence for the
current baseline.

The most valuable immediate test is independent software reproduction from a
clean installation. A future scientific successor may investigate whether
cross-scale persistence can be made more sensitive without reintroducing
false recovery, but it must be separately versioned and evaluated first on
the complete adverse set and then on a new untouched corpus. The historical
9JZO and 4JX2 decisions must remain unchanged.

Exploratory D6c1 should remain deferred under the declared program gate. A
future transition would require an additional committed correct prospective
geometric recovery, zero confident prospective wrong recoveries, and a public
reproducible observation-level reduction. The first downstream experiment
should test only the index bridge and multiplicity-preserving observation
mapping, not intensity prediction.

## 5. Conclusions

A raw-only, cell-free pipeline can construct repeat-certified reciprocal
observations and recover a physical primitive lattice on real public data.
The developed method combines multiscale candidate access with split,
held-out, finite-index, and complexity evidence so that unsupported inputs can
produce explicit abstention rather than a forced lattice.

The current evidence is encouraging but narrow. Baseline B has one correct
prospective recovery, two conservative prospective misses, and no confident
prospective wrong recovery. This establishes real recovery capability, not
broad indexing validity or superiority. The reproducible positive and
negative record defines a concrete boundary for outside testing and future
method development.

## Data and code availability

Source code, frozen method bindings, compact results, evidence tables,
prospective records, and reproduction instructions are available at
<https://github.com/etblink/nfc-crystallography>. The controlling public
software release is `v0.1.1`, commit
`088787ed4927d4c1560e2102d2cb7cd72d073a65`, Git tree
`80411f6268ad987426f376adb961c75215d4437d`. The version-independent Zenodo
concept DOI is <https://doi.org/10.5281/zenodo.21639283>. Public release
`v0.1.1` is the packaging-corrected source of record on GitHub. Its
deterministic source ZIP, wheel, and source-distribution SHA-256 values are,
respectively,
`2f466c1f16d9dd62a5f86cdf8d238f892ec286b43081c746eb4675d4d7e82bee`,
`16d7461c437cf32e3c3fc0c97b29368e37c8f248a79651869083591395213083`,
and
`bafc9b484c7a14f1231ea8ec1d0452aa023dba8bea2bfa494b6feeb3d0153392`.

Raw-data URLs, sizes, and SHA-256 hashes are stored under
`data_manifests/`. The 4JX2 truth-free commitment and post-commitment
diagnostic are preserved under `prospective/baseline_b_transfer_02/`.
The baseline-B archive binding is
`ba18310a04a45c13f1fdf100599c77f2da9fa8ba2f43ec9e942719871b6edf48`.
The controlling 8VTD truth-free decision digest is
`18d8664318f9d27964a4abe179841acf34a2f4bb50e6b836e034ffb154afe249`.
Large raw archives are not redistributed with the software.

## Author contributions

E.T.K. led conceptualization, investigation, methodology development,
software development, validation, evidence curation, and manuscript
preparation.

## Human-AI research disclosure

The research and software were developed through an extended human-AI
collaboration led by E.T.K. AI systems assisted with code generation,
diagnostic design, literature organization, reproducibility packaging, and
drafting. E.T.K. directed the work, set the scientific claim boundaries,
controlled information exposure in prospective tests, verified the released
artifacts, and accepts responsibility for the manuscript. Scientific claims
are based on public inputs, executable code, preserved hashes, and recorded
outcomes rather than on the authority of the conversations that produced
them.

## References

1. Kabsch, W. (2010). XDS. *Acta Crystallographica Section D*, **66**,
   125–132. <https://doi.org/10.1107/S0907444909047337>
2. Winter, G., Waterman, D. G., Parkhurst, J. M., Brewster, A. S.,
   Gildea, R. J., Gerstel, M., Fuentes-Montero, L., Vollmar, M., Michels-Clark,
   T., Young, I. D., Sauter, N. K. & Evans, G. (2018). DIALS: implementation
   and evaluation of a new integration package. *Acta Crystallographica
   Section D*, **74**, 85–97.
   <https://doi.org/10.1107/S2059798317017235>
3. Steller, I., Bolotovsky, R. & Rossmann, M. G. (1997). An algorithm for
   automatic indexing of oscillation images using Fourier analysis.
   *Journal of Applied Crystallography*, **30**, 1036–1040.
   <https://doi.org/10.1107/S0021889897008777>
4. Sauter, N. K., Grosse-Kunstleve, R. W. & Adams, P. D. (2004). Robust
   indexing for automatic data collection. *Journal of Applied
   Crystallography*, **37**, 399–409.
   <https://doi.org/10.1107/S0021889804005874>
5. Gildea, R. J., Waterman, D. G., Parkhurst, J. M., Axford, D., Sutton, G.,
   Stuart, D. I., Sauter, N. K., Evans, G. & Winter, G. (2014). New methods
   for indexing multi-lattice diffraction data. *Acta Crystallographica
   Section D*, **70**, 2652–2666.
   <https://doi.org/10.1107/S1399004714017039>
6. Grosse-Kunstleve, R. W., Sauter, N. K. & Adams, P. D. (2004). Numerically
   stable algorithms for the computation of reduced unit cells. *Acta
   Crystallographica Section A*, **60**, 1–6.
   <https://doi.org/10.1107/S010876730302186X>
7. Brönnimann, C., Eikenberry, E. F., Henrich, B., Horisberger, R., Hülsen,
   G., Pohl, E., Schmitt, B., Schulze-Briese, C., Suzuki, M., Tomizaki, T.,
   Toyokawa, H. & Wagner, A. (2006). The PILATUS 1M detector. *Journal of
   Synchrotron Radiation*, **13**, 120–130.
   <https://doi.org/10.1107/S0909049505038665>
8. Bernstein, H. J. & Hammersley, A. P. (2006). Specification of the
   crystallographic binary file (CBF/imgCIF). In *International Tables for
   Crystallography, Volume G: Definition and exchange of crystallographic
   data*, Chapter 2.3, pp. 37–43.
   <https://doi.org/10.1107/97809553602060000729>
9. Grabowski, M., Cymborowski, M., Porebski, P. J., Osinski, T., Shabalin,
   I. G., Cooper, D. R., Minor, W. & Joachimiak, A. (2016). The Integrated
   Resource for Reproducibility in Macromolecular Crystallography: experiences
   of the first four years. *Acta Crystallographica Section D*, **72**,
   1181–1193. <https://doi.org/10.1107/S2059798316014716>
