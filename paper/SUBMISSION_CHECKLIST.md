# Manuscript submission checklist

This file tracks journal-specific and author-only decisions that should not be
invented in the scientific manuscript.

## Required before submission

- [ ] Select the target journal and apply its article template, word limit,
  reference style, and data-policy requirements.
- [ ] Confirm the corresponding-author email and whether an ORCID should be
  listed.
- [ ] Confirm the final affiliation wording.
- [ ] Add the author's funding statement.
- [ ] Add the author's competing-interests declaration.
- [ ] Confirm whether the target journal requires a separate data-availability,
  code-availability, ethics, or AI-use form.
- [ ] Confirm that the Zenodo `v0.1.1` files are published under concept DOI
  `10.5281/zenodo.21639283`; until then, cite the exact GitHub release as the
  `v0.1.1` source of record.
- [ ] Record any Software Heritage snapshot identifier after archival
  completes.
- [ ] Decide whether the detailed evidence table remains in the main text or
  moves to supplementary material.
- [ ] Add figures only when they communicate a relationship not already clear
  from the tables. At minimum, consider a pipeline schematic and the D4.5
  causal-ablation comparison.
- [ ] Ask one crystallographer or diffraction-software developer to review the
  method definitions and claim boundary.

## Scientific checks

- [ ] Keep baseline B unchanged and controlling.
- [ ] Keep candidate C labeled as an experimental compatible alternative with
  no demonstrated incremental value.
- [ ] Preserve 8VTD as the single correct prospective recovery.
- [ ] Preserve 9JZO and 4JX2 as conservative missed recoveries; do not count
  latent or single-scale physical families as recovered decisions.
- [ ] Preserve 4G2A, 6CKT, and 6TPI as pre-execution exclusions with no NFC
  result.
- [ ] State that three prospective corpora do not support a general success-rate
  estimate.
- [ ] State that independent external reproduction and D6c1-D7 transfer are
  not established.
- [ ] Keep the crystallographic claim independent of NFC cosmology or
  Theory-of-Everything claims.

## Release validation

- [ ] Run `ruff check .`.
- [ ] Run `pytest`.
- [ ] Run `nfc-cryst verify-methods`.
- [ ] Run the compact 8VTD replay and confirm zero nonfloating mismatches.
- [ ] Build wheel and source distribution.
- [ ] Install the wheel in a clean Python 3.10 environment and a clean Python
  3.12 environment.
- [ ] Verify the raw-builder imports from the installed wheel.
- [ ] Record the final manuscript commit and Git tree in the submission
  materials.
