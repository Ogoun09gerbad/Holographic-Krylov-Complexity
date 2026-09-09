# Connecting GitHub → Zenodo and minting the DOI

This repository does not yet have a GitHub URL or a Zenodo DOI — both require manual steps in your own accounts that cannot be performed on your behalf. Follow these steps after the repository has been pushed to GitHub.

## 1. Create the GitHub repository

```bash
gh repo create Holographic-Krylov-Complexity --public --source=. --remote=origin
# or manually on github.com, then:
git remote add origin https://github.com/<your-username>/Holographic-Krylov-Complexity.git
git push -u origin main
```

## 2. Replace the placeholders

Search the repository for `REPLACE_WITH_GITHUB_USERNAME` and replace it with your actual GitHub username/organization in:
- `README.md`
- `CITATION.cff`
- `report/Holographic_Krylov_Complexity_MIORPA_Report.tex` (the "Code and reproducibility" paragraph near the end of the introduction)

Commit this change before tagging the release.

## 3. Tag and publish the v1.0.0 release

```bash
git tag -a v1.0.0 -m "v1.0.0 — version corresponding to the MIORPA 2026 submitted report"
git push origin v1.0.0
```

Then, on GitHub: **Releases → Draft a new release → choose tag v1.0.0** → title it `v1.0.0` → publish.

## 4. Connect the repository to Zenodo (one-time setup)

1. Go to <https://zenodo.org>, sign in with your GitHub account.
2. Go to **GitHub** in your Zenodo account settings, find `Holographic-Krylov-Complexity` in the repository list, and toggle it **on**. (If it doesn't appear, click "Sync now".)
3. Zenodo will now automatically archive every future GitHub release of this repository.

## 5. Trigger the archive for v1.0.0

If you enabled the Zenodo webhook *before* creating the v1.0.0 release, it is archived automatically. If you enabled it *after*, either:
- delete and re-publish the v1.0.0 GitHub release, or
- create a v1.0.1 release once the webhook is active, and treat that as the first archived version.

## 6. Retrieve the DOI

Once archived, Zenodo assigns a DOI (e.g. `10.5281/zenodo.XXXXXXX`) shown on the deposit's Zenodo page and on a badge you can copy from **Zenodo → your deposit → "DOI" badge markdown**.

## 7. Insert the DOI everywhere it is currently a placeholder

Replace `REPLACE_WITH_ZENODO_DOI` / `TO_BE_REPLACED_AFTER_ZENODO_PUBLICATION` with the real DOI in:
- `README.md` (badge at the top, BibTeX entry, "Links" section)
- `CITATION.cff` (`doi:` field)
- `report/Holographic_Krylov_Complexity_MIORPA_Report.tex` (the `\url{https://doi.org/...}` in the "Code and reproducibility" paragraph)

Recompile the report PDF after this change and commit both the `.tex` and the updated PDF.

## Workflow summary

```text
GitHub repository  →  GitHub release v1.0.0  →  Zenodo archive  →  DOI  →  DOI inserted into CITATION.cff and the article
```

No DOI has been fabricated anywhere in this repository — every occurrence is the literal placeholder string `REPLACE_WITH_ZENODO_DOI` or `TO_BE_REPLACED_AFTER_ZENODO_PUBLICATION`, so a simple repository-wide search will find every place that needs the real value.
