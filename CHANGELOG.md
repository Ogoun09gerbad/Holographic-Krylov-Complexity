# Changelog

## [1.0.0] — 2026-09-09

First public release, corresponding exactly to the version of the report submitted for MIORPA 2026.

### Included
- Final report source, bibliography, and style file (`report/Holographic_Krylov_Complexity_MIORPA_Report.tex`, `Krylov.bib`, `JHEP.bst`) and a freshly recompiled PDF.
- All three verified computational components: `ads3/`, `quiver_comparison/`, `blackstring_quiver/`, each with their code, generated figures, and CSV/table outputs re-run and confirmed against the report on this date.
- `notebooks/Notebook_Krylov_Final.ipynb`, the exploratory notebook predating the `quiver_comparison/` scripts.
- `report/supplementary/Notes.tex`, working derivation notes cited by name in code comments.
- Full documentation set: `README.md`, `REPRODUCIBILITY.md`, `REPOSITORY_AUDIT.md`, `SCIENTIFIC_ISSUES.md`, `CITATION.cff`, `LICENSE`, `ZENODO_RELEASE.md`, sub-directory `README.md` files.

### Changed relative to the original project directory (build/organization only — no science altered)
- `report/Holographic_Krylov_Complexity_MIORPA_Report.tex`: updated `\graphicspath` and eight `\IfFileExists{...}` figure paths from `quiver_comparison/figures/...` / `blackstring_quiver/figures/...` / bare filenames to `../quiver_comparison/figures/...` / `../blackstring_quiver/figures/...` / `../ads3/...` / `../figures/...`, to match this repository's layout (report now lives in `report/`, one level below the code directories). No equation, numerical value, or scientific claim was changed.
- Added a "Code and reproducibility" paragraph at the end of the introduction, with placeholder GitHub/Zenodo links to be filled in after publication (see `REPOSITORY_AUDIT.md`, ACTION REQUIRED).
- Recompiled the PDF from the relocated source (`pdflatex` ×3 + `bibtex`); zero undefined references in the final pass.

### Excluded from this repository (left untouched in the original project directory)
- Presentation materials: `krylov_beamer*.tex/pdf`, `QA_prep_soutenance.tex`, `fiche_1page_reunion.tex`, `speech_notes_*.md`, `oxford_assets/` — out of scope for a code-reproducibility repository.
- Earlier report drafts `Report.tex`, `Plan.tex`, and the apparently-unused `Krylov_MIORPA.bib` — see `REPOSITORY_AUDIT.md` for why.
- Two stale figure files in `blackstring_quiver/figures/` not reproduced by the current script (`02_J_independent_equilibrium.png`, `03_barrier_crossing_trajectories.png`).
- Local virtual environments `krylov_env/`, `.venv/` — replaced by `requirements.txt`.

### Known gaps (see `REPOSITORY_AUDIT.md` and `SCIENTIFIC_ISSUES.md`)
- `figures/ct_example.png` (Fig. `fig:toy-ct`) has no recovered generating script.
- GitHub repository creation and Zenodo archival are manual follow-up steps, not yet performed.
