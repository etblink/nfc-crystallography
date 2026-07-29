# Manuscript submission checklist

This file tracks journal-specific and author-only decisions that should not be
invented in the scientific manuscript.

## Required before submission

- [x] Select *Journal of Applied Crystallography*, Research Paper, as the
  target. Record its scope, category choice, length guidance, reference style,
  data policy, and publication route in `paper/JOURNAL_PLAN.md`.
- [ ] Convert the Markdown manuscript to the current IUCr Word or LaTeX
  template.
- [x] Add the corresponding-author email and ORCID.
- [x] Confirm the affiliation as “Independent researcher (solo,
  AI-assisted), Las Vegas, Nevada, USA”.
- [x] Add the author's self-funding statement.
- [x] Add the author's competing-interests declaration.
- [x] Include explicit data-and-code availability and human–AI-use
  disclosures. Recheck the live submission form for any corresponding
  metadata fields when submitting.
- [x] Draft the journal covering letter.
- [ ] Confirm that the Zenodo `v0.1.1` files are published under concept DOI
  `10.5281/zenodo.21639283`; until then, cite the exact GitHub release as the
  `v0.1.1` source of record.
- [ ] Record any Software Heritage snapshot identifier after archival
  completes.
- [ ] Decide whether the detailed evidence table remains in the main text or
  moves to supplementary material.
- [ ] Add the graphical abstract or thumbnail requested by the journal. Add
  other figures only when they communicate a relationship not already clear
  from the tables; at minimum, consider a pipeline schematic and the D4.5
  causal-ablation comparison.
- [ ] Ask one crystallographer or diffraction-software developer to review the
  method definitions and claim boundary.
- [ ] Confirm immediately before submission that the manuscript is original
  and is not under consideration elsewhere.

## Scientific checks

- [x] Keep baseline B unchanged and controlling.
- [x] Keep candidate C labeled as an experimental compatible alternative with
  no demonstrated incremental value.
- [x] Preserve 8VTD as the single correct prospective recovery.
- [x] Preserve 9JZO and 4JX2 as conservative missed recoveries; do not count
  latent or single-scale physical families as recovered decisions.
- [x] Preserve 4G2A, 6CKT, and 6TPI as pre-execution exclusions with no NFC
  result.
- [x] State that three prospective corpora do not support a general success-rate
  estimate.
- [x] State that independent external reproduction and D6c1-D7 transfer are
  not established.
- [x] Keep the crystallographic claim independent of NFC cosmology or
  Theory-of-Everything claims.

## Release validation

- [x] Run `ruff check .`.
- [x] Run `pytest` (14/14 passed).
- [x] Run `nfc-cryst verify-methods` (21/21 source bindings passed).
- [x] Run the compact 8VTD replay and preserve
  `LATTICE_RECOVERED / ONE_PERSISTENT_FAMILY_SURVIVES` with decision digest
  `18d8664318f9d27964a4abe179841acf34a2f4bb50e6b836e034ffb154afe249`.
- [x] Build wheel and source distribution.
- [x] Install and verify the wheel under Python 3.10 and Python 3.12 in GitHub
  Actions run `30415938565`.
- [x] Verify the raw-builder imports and all 21 frozen method sources from an
  isolated installed wheel under Python 3.12.
- [ ] Record the final manuscript commit and Git tree in the submission
  materials.
