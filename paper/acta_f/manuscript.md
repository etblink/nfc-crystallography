# Raw-image lattice recovery with repeat-certified spots and explicit abstention

**Evan Thomas Kotler**

Independent researcher (solo, AI-assisted), Las Vegas, Nevada, USA

Correspondence: evantkotler@gmail.com

ORCID: <https://orcid.org/0009-0004-5840-4443>

**Article type:** Methods Communication

**Synopsis:** A cell-blind method groups repeated diffraction-spot detections,
searches for primitive reciprocal lattices at several length scales, and
reports either a supported lattice or an explicit abstention. Tests on public
macromolecular diffraction images show real recovery capability together with
clear sensitivity limits.

## Abstract

Indexing rotation diffraction images assigns spot positions to a reciprocal
lattice and determines its orientation. An exploratory method was developed
to make this decision from raw images without using a supplied unit cell,
orientation matrix, symmetry assignment or processed reflection file.
Frame-local spot detections are associated across adjacent images and
consolidated into three-dimensional reciprocal-space centroids. Primitive
lattice candidates are then generated at three reciprocal-length scales and
must survive independent half-scan fits, cross-half prediction, shifted-grid
controls and tests of related coarser or finer lattices. The method can return
`LATTICE_RECOVERED`, `AMBIGUOUS_LATTICE` or `INSUFFICIENT_SIGNAL`.

On public 9JQ9 images, 58,528 detections matched cell-blind DIALS spots within
one pixel and their integrated intensities had a Spearman correlation of
0.903. Controlled 6GN3 tests showed that consolidating repeat detections into
reciprocal centroids, rather than filtering alone, caused recovery by the
unchanged lattice search. A fixed development census returned all 32 required
recovery or abstention outcomes with no incorrect recovery. In three
prospective public-data tests, the released method correctly recovered the
primitive lattice for 8VTD and abstained for 9JZO and 4JX2. Post-decision
scoring found a physical primitive candidate in each missed case, but neither
was promoted retrospectively. The prospective record is therefore one correct
recovery, two conservative misses and no confident incorrect recovery.

These results establish a raw-derived lattice-candidate and abstention
framework with genuine but limited transfer. They do not establish superiority
to established indexing software, a production-ready indexer or transfer to
reflection intensities and structure determination.

**Keywords:** diffraction indexing; reciprocal lattice; rotation data; spot
finding; macromolecular crystallography; abstention; reproducible research

## 1. Introduction

Before intensities from a rotation diffraction experiment can be merged or
used for structure determination, the diffraction spots must be assigned
integer indices. Indexing finds a reciprocal lattice and an orientation that
explain the measured spot positions. Mature packages such as XDS and DIALS
perform this task within comprehensive data-processing systems (Kabsch, 2010;
Winter et al., 2018). Established indexing strategies include Fourier analysis
of reciprocal-space difference vectors, real-space grid searches and methods
for multiple lattices (Steller et al., 1997; Sauter et al., 2004; Gildea et
al., 2014).

The present study asks a deliberately narrower question. Can a representation
constructed directly from raw diffraction images produce a reproducible
primitive-lattice decision without being given a conventional indexing
solution? Here, *primitive lattice* means the smallest repeating reciprocal
grid supported by the observations. A method should not be judged only by
whether it can fit a known positive example. It should also avoid selecting a
dense but physically unsupported grid and should decline to report a lattice
when the data do not distinguish one explanation.

This requirement arose from an early version of the method. That version could
recover clean synthetic lattices, but it also assigned approximately 96%
formal support to a deliberately nonlattice point cloud. High support was
therefore not sufficient evidence of a physical lattice. Subsequent public
real-data and synthetic experiments isolated two additional problems. First,
one physical diffraction spot can be detected on several adjacent rotation
images, so frame-local detections need to be consolidated before lattice
search. Second, fixed reciprocal-length limits can hide modes associated with
large direct-space unit cells.

The resulting public method is referred to as **baseline B**. It has three
implementation stages (Fig. 1). D4 detects spot-like signal independently on
each image. D4.5 links compatible detections across adjacent images and
replaces each accepted repeat path with a reciprocal-space centroid. D5
generates primitive-lattice candidates at several scales and applies
independent-data and alternative-model checks. It reports a lattice only when
one persistent family survives; otherwise it reports ambiguity or insufficient
signal.

The aim of this Methods Communication is to explain that workflow in terms
accessible to structural biologists and to report its positive and negative
evidence together. The method is not presented as a replacement for DIALS,
XDS or other production indexing programs. Its narrower contribution is a
raw-derived candidate-and-abstention framework whose complete prospective
record currently contains one correct recovery and two missed recoveries.

![Baseline B starts from native rotation images. Frame-local spot detections
(D4) are linked across adjacent images and consolidated into reciprocal-space
centroids (D4.5). A multiscale search generates primitive-lattice candidates
(D5). Recovery is authorized only when split, held-out, finite-index and
complexity checks leave one persistent family; otherwise the method
abstains.](../figures/figure_1_pipeline.png){width=88mm}

## 2. Materials and methods

### 2.1 Public images and evidence roles

The reported experiments used public macromolecular rotation datasets from
the Integrated Resource for Reproducibility in Macromolecular Crystallography
and associated ProteinDiffraction repositories (Grabowski et al., 2016).
The public raw-image path accepts native PILATUS 1.2 MiniCBF data with
`x-CBF_BYTE_OFFSET` compression. Detector dimensions, pixel size, distance,
beam centre, axes, wavelength and image angle are taken from the image model.
Images are not resampled or cropped, and pixel values are not scaled.

Datasets were assigned roles according to information exposure. Development
datasets could be inspected and used to change the method. They therefore
provide mechanistic or retrospective evidence, not out-of-sample validation.
For a prospective test, the corpus and fixed software identity were recorded
before execution. Conventional indexing established that the full scan and
independent scan halves had stable solutions, but the conventional unit cell
and orientation were concealed from the experimental method. The method's
decision was recorded and hash-bound before the conventional solution was
revealed for scoring. A dataset excluded before method execution did not count
as a method outcome.

### 2.2 Frame-local spot detection

For each image, D4 estimates background and robust scale separately in
32 × 32-pixel tiles using the median, median absolute deviation and a
Poisson-like lower scale bound. Fixed high-threshold pixels seed connected
components, which grow through lower-threshold neighbouring pixels. Very small
or weak components are rejected. Accepted components are represented by
signal-weighted pixel centroids and by deterministic pixel and reciprocal-space
enclosures.

This stage uses neither a unit cell nor a crystal orientation. The enclosures
describe the spatial extent allowed by the construction and are used only for
association; they are not calibrated confidence intervals. Detector rays are
mapped to reciprocal coordinates using the native beam, detector and scan
geometry. The implementation uses reciprocal units of cycles per ångström.

### 2.3 Repeat-certified reciprocal spots

During rotation, the same reciprocal spot may leave signal on consecutive
images. Treating every frame-local detection as an independent
three-dimensional point can fragment or duplicate the evidence supplied to a
lattice search. D4.5 therefore processes each acquisition series separately
and follows the sign of its rotation increment. Associations are allowed only
between consecutive images whose reciprocal enclosures overlap. Ambiguous
many-to-one associations are rejected rather than resolved using a lattice
model.

Each accepted path is replaced by a signal-weighted reciprocal centroid. The
centroid must remain within every contributing enclosure, and the association
must have a unique contiguous interpretation. Full-scan, first-half and
second-half feeds are constructed independently. No conventional cell,
orientation, symmetry or processed reflection is used. The result is a
deterministically ordered set of repeat-certified three-dimensional
observations.

### 2.4 Multiscale candidate search

D5 searches difference vectors between reciprocal centroids because
translations between lattice points contain information about the lattice
periodicity. The fixed candidate generator balances observations across
acquisition series, collects local difference vectors, bins them at
0.003 Å⁻¹ and forms independent triples from recurring difference-vector
modes. Candidate bases are refined by assigning nearby reciprocal points to
integer coordinates and refitting.

The search is repeated with lower difference-vector limits of 0.012, 0.009 and
0.006 Å⁻¹. These values are four, three and two times the fixed bin width.
Each scale is fitted independently to the full scan and to both scan halves.
A family must recur at two or more scales before it can be reported as
recovered. A strong candidate found at only one scale can remain visible in
the diagnostic record but cannot authorize recovery.

### 2.5 Evidence tests and abstention

Four tests separate candidate generation from the final decision.

1. **Split reproducibility.** The full-scan and two half-scan fits must describe
   the same primitive lattice after allowing integer changes of lattice basis.
2. **Held-out prediction.** A basis fitted to one half must predict periodic
   positions in the other half. Its zero-origin grid is compared with 63
   equally dense grids shifted by quarter-cell fractions.
3. **Related-lattice alternatives.** Seven index-two primitive completions are
   tested so that a convincing finer primitive lattice is not hidden by an
   index-two alias.
4. **Complexity control.** Related coarser lattices of integer index 2–8 are
   tested. If a lower-complexity model retains nearly all of the selected
   support, recovery is withheld.

These checks answer different questions. Split disagreement detects
scan-dependent solutions; shifted-grid controls test whether the fitted origin
is specific; the index-two audit exposes competing primitive explanations; and
the complexity test limits the tendency of dense reciprocal grids to collect
points by chance. Fixed precedence maps the evidence to one of three outcomes:
`LATTICE_RECOVERED`, `AMBIGUOUS_LATTICE` or `INSUFFICIENT_SIGNAL`.

An experimental comparator, candidate C, replaces only the hard
origin-neighbourhood count with a continuous Cauchy affinity of fixed
0.003 Å⁻¹ scale. It preserved all development decisions and tied baseline B
exactly in one prospective comparison. Because it supplied no incremental
evidence, it remains experimental and is not the reported method.

### 2.6 Conventional scoring and reproducibility

Prospective candidates were scored only after decision commitment. Conventional
full- and half-scan indexing had to assign at least half of the accepted spots
and agree in reciprocal metric and orientation. The orientation comparison
minimizes disagreement over a finite set of proper unimodular basis changes;
this avoids rejecting equivalent indexing representations merely because they
use different bases. Experimental and conventional bases were then compared
for primitive equivalence by the nearest integer basis transformation.

The Python implementation, fixed method files, tests, compact result records,
prospective commitments and raw-data manifests are public. Release `v0.1.1`
is a packaging-only correction that restored two frozen support modules
omitted from `v0.1.0`; no scientific rule or result changed. A compact 8VTD
replay reproduces the archived decision with no nonfloating mismatch and a
maximum floating-point difference of 2.83 × 10⁻¹⁵.

## 3. Results

### 3.1 The upstream stages capture diffraction information

The 9JQ9 dataset provided a detector and specimen distinct from the initial
development case. D4 produced 58,528 mutual pairs within one pixel of
cell-blind DIALS spot centroids. Their integrated intensities had Spearman
correlation 0.903. This shows that D4 responds to real diffraction signal. It
does not show that D4 is more accurate or informative than DIALS.

The role of repeat consolidation was tested on two 6GN3 sweeps while the
downstream lattice search was held unchanged. Three inputs were compared:
unaggregated frame-local detections, only the original members of
repeat-certified paths, and the consolidated reciprocal centroids. The first
two inputs failed on both sweeps. Arithmetic and signal-weighted centroids
recovered the physical primitive lattice on both full sweeps and independently
constructed halves. Repeat filtering alone was therefore insufficient;
moving repeated detections to a shared three-dimensional centroid caused the
observed repair.

A separate scan-direction test found that an ascending-only traversal produced
no repeat objects for the negative-increment 9Z6F scans. The released
scan-local, increment-aware traversal produced 17,880 repeat objects and
recovered the historical primitive lattice for the full dataset and both
scan-configuration partitions.

### 3.2 Development controls test recovery and refusal

The predecessor D5 objective recovered four idealized synthetic lattices but
failed all six required adverse abstentions. Most strikingly, a pure
nonlattice point cloud received 95.625% formal support. This result motivated
the multiscale persistence, related-lattice and complexity checks.

The developed baseline produced the required result for all 32 cases in its
fixed development census: 16 recoveries, two ambiguity decisions and 14
insufficient-signal decisions, with no incorrect recovery. The set included
synthetic primitive lattices, missing-data and masking perturbations,
competing lattices, sparse superlattice signal and nonlattice backgrounds, as
well as public real-data controls.

The real-data lineage recovered primitive lattices on 9Z6F, 6GN2 and two 6GN3
sweeps. It also recovered 6MFU after that corpus entered development. The 6MFU
diagnosis showed why both multiscale access and the index-two completion audit
were needed: a required reciprocal vector lay below the original
0.012 Å⁻¹ limit, and the physical primitive completion was in the direction
not examined by a coarsening-only test.

Not every public dataset became a recovery after these repairs. 9JQ9 remained
split-unstable, and its public raw-to-deposited provenance could not be fully
resolved. 5V0G and the five-series 6W61 dataset also remained insufficient.
These failures helped define sensitivity limits but were not converted into
positive results.

### 3.3 Prospective public-data outcomes

Table 1 gives the complete prospective record for the released method after
development. A plausible latent candidate does not count as recovery unless
the committed decision was `LATTICE_RECOVERED`.

**Table 1. Prospective real-data record.**

| Dataset | Committed decision | Post-decision assessment |
| --- | --- | --- |
| 8VTD | `LATTICE_RECOVERED` | Correct physical primitive lattice |
| 9JZO | `AMBIGUOUS_LATTICE` | Conservative miss; a latent survivor family was primitive-equivalent |
| 4JX2 | `INSUFFICIENT_SIGNAL` (`NO_PERSISTENT_FAMILY`) | Conservative miss; a primitive-equivalent family occurred at one scale |

The 8VTD result is the present method's strongest transfer evidence. The
method and inputs were fixed before the conventional result was consulted, and
the committed family was primitive-equivalent to the conventional full- and
split-scan solutions.

For 9JZO, baseline B and experimental candidate C returned the same
`AMBIGUOUS_LATTICE` decision and selected the same latent families.
Post-decision scoring found their common surviving family to be
primitive-equivalent. The committed result remains a missed recovery, and the
tie provides no reason to replace baseline B with C.

For 4JX2, 540 of 542 native images entered the fixed feeds. D4.5 produced
31,668 full-scan, 16,024 first-half and 15,647 second-half reciprocal
centroids. Baseline B returned `INSUFFICIENT_SIGNAL` because no family
persisted at two scales. After commitment, the direct family at
0.006 Å⁻¹ was found to be primitive-equivalent to all three conventional
solutions. Because it was present at only one tested scale, the abstention was
not changed.

The prospective record is therefore one correct recovery, two conservative
misses and no confident incorrect recovery. Three datasets are too few to
estimate a general success or error rate. The record instead shows that real
out-of-sample recovery is possible and that the current persistence rule can
also reject a physically correct candidate.

## 4. Discussion

The main result is a separation between *finding a plausible lattice* and
*earning permission to report it*. Repeat-certified reciprocal centroids and
multiscale difference-vector modes contain genuine lattice information: they
recover several public development lattices and one untouched prospective
lattice. At the same time, the nonlattice and competing-lattice controls show
why raw support cannot safely be used alone.

The layered decision is intentionally conservative. Split reproducibility,
held-out prediction, related-lattice alternatives and complexity control
reject distinct failure modes. That conservatism prevented a confident wrong
answer in the three prospective tests, but it also produced two missed
recoveries. In 9JZO the physical family survived but remained ambiguous; in
4JX2 it appeared at only one scale. These outcomes identify candidate
persistence and final discrimination as the current sensitivity boundary.
They are not hidden successes.

Several limitations constrain the claim. The public raw path presently covers
native PILATUS 1.2 MiniCBF rather than the detector breadth of established
suites. The development and prospective collections are small and not
representative samples of macromolecular diffraction experiments. The study
does not provide a broad speed or accuracy benchmark against DIALS, XDS or
other indexers. The reciprocal enclosures used for repeat association are
deterministic construction bounds, not calibrated measurement uncertainties.
Candidate generation and the related-lattice audits are finite and do not
enumerate every lattice. Large raw archives remain external public
dependencies. Independent clean-install reproduction by an outside group has
not yet been reported.

Most importantly, this work stops at lattice geometry. It does not establish
observation-level index transfer, preserved multiplicity, intensity analysis,
symmetry assignment, structure solution or refinement. A future downstream
study should begin only after further prospective geometric recovery and
should first test the narrow bridge between raw-derived lattice indices and a
reproducible conventional observation-level reduction.

The appropriate current conclusion is therefore limited but positive. A
cell-blind, raw-image representation can support primitive-lattice recovery
with explicit abstention, and that ability has transferred prospectively once.
The same public record also documents two conservative misses. This balance
makes the method falsifiable and gives outside users a concrete software and
data record to reproduce, criticize and extend.

## 5. Data and code availability

Source code, frozen method bindings, tests, evidence tables, prospective
records and raw-data manifests are available at
<https://github.com/etblink/nfc-crystallography>. The controlling public
software release is `v0.1.1`. Its version DOI is
<https://doi.org/10.5281/zenodo.21654414>, linked to the all-versions concept
DOI <https://doi.org/10.5281/zenodo.21639283>. Large raw archives are not
redistributed; public locations, byte sizes and SHA-256 hashes are recorded in
the repository.

## Acknowledgements

The public raw diffraction repositories and the developers of DIALS, dxtbx,
XDS and the broader open crystallographic software ecosystem made this
evaluation possible.

## Author contributions

E.T.K. led conceptualization, investigation, methodology and software
development, validation, evidence curation and manuscript preparation.

## Human–AI research disclosure

The research and software were developed through an extended human–AI
collaboration led by E.T.K. AI systems assisted with code generation,
diagnostic design, literature organization, reproducibility packaging and
drafting. E.T.K. directed the work, controlled information exposure in
prospective tests, verified the released artifacts, checked the manuscript and
accepts responsibility for its content. AI systems are not authors.

## Funding information

This work was self-funded by the author and received no external funding.

## Conflict of interest

The author declares no competing interests.

## References

Gildea, R. J., Waterman, D. G., Parkhurst, J. M., Axford, D., Sutton, G.,
Stuart, D. I., Sauter, N. K., Evans, G. & Winter, G. (2014). New methods for
indexing multi-lattice diffraction data. *Acta Crystallographica Section D*,
**70**, 2652–2666. <https://doi.org/10.1107/S1399004714017039>

Grabowski, M., Cymborowski, M., Porebski, P. J., Osinski, T., Shabalin, I. G.,
Cooper, D. R., Minor, W. & Joachimiak, A. (2016). The Integrated Resource for
Reproducibility in Macromolecular Crystallography: experiences of the first
four years. *Acta Crystallographica Section D*, **72**, 1181–1193.
<https://doi.org/10.1107/S2059798316014716>

Kabsch, W. (2010). XDS. *Acta Crystallographica Section D*, **66**, 125–132.
<https://doi.org/10.1107/S0907444909047337>

Sauter, N. K., Grosse-Kunstleve, R. W. & Adams, P. D. (2004). Robust indexing
for automatic data collection. *Journal of Applied Crystallography*, **37**,
399–409. <https://doi.org/10.1107/S0021889804005874>

Steller, I., Bolotovsky, R. & Rossmann, M. G. (1997). An algorithm for
automatic indexing of oscillation images using Fourier analysis. *Journal of
Applied Crystallography*, **30**, 1036–1040.
<https://doi.org/10.1107/S0021889897008777>

Winter, G., Waterman, D. G., Parkhurst, J. M., Brewster, A. S., Gildea, R. J.,
Gerstel, M., Fuentes-Montero, L., Vollmar, M., Michels-Clark, T., Young,
I. D., Sauter, N. K. & Evans, G. (2018). DIALS: implementation and evaluation
of a new integration package. *Acta Crystallographica Section D*, **74**,
85–97. <https://doi.org/10.1107/S2059798317017235>
