# Manuscript submission checklist

This file tracks journal-specific and author-only decisions that should not be
invented in the scientific manuscript.

## Required before submission

- [x] Select *Journal of Applied Crystallography*, Research Paper, as the
  target. Record its scope, category choice, length guidance, reference style,
  data policy, and publication route in `paper/JOURNAL_PLAN.md`.
- [x] Prepare a single-column, double-spaced Word manuscript consistent with
  the current IUCr preprint requirements. The optional hosted DOCX template
  was unavailable behind an interactive verification page, so the distributable
  Word source was built directly from the published formatting requirements.
- [x] Add the corresponding-author email and ORCID.
- [x] Confirm the affiliation as “Independent researcher (solo,
  AI-assisted), Las Vegas, Nevada, USA”.
- [x] Add the author's self-funding statement.
- [x] Add the author's competing-interests declaration.
- [x] Include explicit data-and-code availability and human–AI-use
  disclosures. Recheck the live submission form for any corresponding
  metadata fields when submitting.
- [x] Draft the journal covering letter.
- [x] Confirm the four exact Zenodo `v0.1.1` files under version DOI
  `10.5281/zenodo.21654414`, linked to concept DOI
  `10.5281/zenodo.21639283`.
- [ ] Record the Software Heritage snapshot identifier when save request
  `2402658` completes. This is archival redundancy and is not a submission
  blocker because the exact release is already public through GitHub and
  Zenodo.
- [x] Keep both evidence tables in the main text. Table 1 defines the
  development/prospective/exclusion roles needed to interpret every result,
  and Table 2 makes the complete prospective record visible. No supplementary
  file is warranted solely to duplicate repository JSON records.
- [x] Add a 600-dpi pipeline figure to the manuscript and a separate
  600-dpi graphical abstract.
- [x] Publicly request an independent clean-install/replay review through
  GitHub issue #6 and external posts. No response has yet been received;
  independent reproduction remains explicitly unestablished and is not being
  treated as a submission blocker.
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
- [ ] Record the merge commit and Git tree after the manuscript PR is merged;
  then update the private submission log or journal form without changing the
  submitted scientific content.
