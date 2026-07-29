# Journal of Applied Crystallography submission package

This directory is generated from the repository manuscript and contains the
ordinary files needed for a *Journal of Applied Crystallography* Research Paper
submission.

## Files to upload

1. `JAC_submission_manuscript.docx` — single-column, double-spaced Word
   manuscript with the pipeline figure embedded.
2. `figure_1_pipeline.png` — separate 2400 × 1350, 600-dpi figure file.
3. `graphical_abstract.png` — separate 2400 × 1350, 600-dpi graphical
   abstract/thumbnail.
4. `JAC_cover_letter.docx` — covering letter.

No supporting-information file is supplied. Both compact evidence tables are
central to the interpretation and remain in the article; detailed
machine-readable evidence is already public in the repository and Zenodo
release.

`SHA256SUMS.txt` binds the finalized handoff files. It is for local
verification and is not itself a journal upload.

## Files for form entry and handoff

- `abstract.txt`
- `synopsis.txt`
- `keywords.txt`
- `suggested_referees.md`
- `submission_metadata.json`
- `SUBMISSION_HANDOFF.md`
- `SHA256SUMS.txt`

## Rebuild

From the repository root:

```bash
python paper/tools/build_submission_package.py
```

The generated Word files should be rendered and visually inspected after any
source edit. The build script does not change scientific code or method
definitions.
