# Acta Crystallographica Section F rewrite

This directory contains the general-reader rewrite prepared after the
*Journal of Applied Crystallography* scope decision and the invitation to seek
informal advice from the Managing Editor of *Acta Crystallographica Section F*.

The rewrite is a **Methods Communication**. It preserves the released
scientific method and evidence record:

- baseline B remains controlling and unchanged;
- candidate C remains an experimental compatible alternative with no
  demonstrated incremental value;
- 8VTD remains the single correct prospective recovery;
- 9JZO and 4JX2 remain conservative missed recoveries;
- no D6c1–D7, production-indexer, independent-validation or broader physical
  claim is added.

Build the Word manuscript from the repository root with:

```bash
"$CODEX_PRIMARY_RUNTIME_PYTHON" \
  paper/tools/build_acta_f_methods_communication.py
```

The journal-facing files are written under `paper/acta_f/submission/`.
