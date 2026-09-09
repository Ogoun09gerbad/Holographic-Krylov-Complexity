# Repository audit

Honest audit performed 2026-09-09, from a clean copy of the source project, before publishing this repository. Every script listed under PASS was actually executed (not just read) in an isolated temporary directory using `C:\Python314\python.exe` (numpy 2.5.2, scipy 1.18.0, sympy 1.14.0, matplotlib 3.11.1) and its output compared to the numbers/figures quoted in `Holographic_Krylov_Complexity_MIORPA_Report.tex`.

## PASS

- `ads3/ads3_validation.py` runs end-to-end; all 6 symbolic checks PASS; worst numeric relative error $1.04\times10^{-7}$, matches the report's "at most $1.0\times10^{-7}$".
- `quiver_comparison/run_quiver_comparison.py` runs end-to-end; naive intercepts reproduced exactly as $0.0877, 0.2630, 0.2307, 0.6920$ (report Table, Sec. `sec:quivers`); $N$-scaling self-check passes to machine precision; exact identity `max|Pr-Ht|` at the $10^{-12}$ level.
- `quiver_comparison/run_intercept_plot.py` runs end-to-end, produces `09_naive_minus_routhian_intercept.png`.
- `quiver_comparison/run_parameter_family_scan.py` runs end-to-end, produces `09_f1_shape_vs_P.png` and `10_scan_vs_P.png` and `parameter_family_scan.csv`.
- `quiver_comparison/sympy_verify_quiver2.py` runs end-to-end; symbolic re-derivation matches `geometry.py` exactly (difference `= 0` for all three quiver-2 pieces).
- `blackstring_quiver/run_blackstring_quiver.py` runs end-to-end (imports `../quiver_comparison/routhian.py` correctly by relative path); Table A–D reproduced; $\Delta z\approx4.49$–$4.52\times10^{-7}$, matches "order $10^{-7}$" claim; `max|F1/F3-4| = 0.000e+00`.
- `notebooks/Notebook_Krylov_Final.ipynb` executes end-to-end with no errors via `jupyter nbconvert --execute` (300s timeout, not hit).
- `report/Holographic_Krylov_Complexity_MIORPA_Report.tex` compiles cleanly from its new `report/` location (`pdflatex` ×3 + `bibtex`, TinyTeX) with **zero** undefined references/citations after the final pass. Output: 33-page PDF.
- No hardcoded machine-specific absolute paths found in any `.py` file in the repository (checked with a recursive grep for `C:\Users`, `/c/Users`, `C:/Users`); all scripts resolve paths via `os.path.dirname(os.path.abspath(__file__))`.
- No randomness is used anywhere in the numerical pipeline (all integrations are deterministic ODE solves), so no seeding was required or added.
- No LICENSE or README existed anywhere in the original project directory tree (verified by recursive search) before this audit — MIT was chosen and applied fresh, not silently overriding an existing choice.

## WARNING

- `figures/ct_example.png` (Fig. `fig:toy-ct`, the $3\times3$ toy-model figure in Sec. 2) has **no surviving generating script** anywhere in the original project directory — not in any `.py` file, not in the notebook (which contains no `savefig` calls at all). The underlying numerical claim in the text ($a_0=2,\ b_1=1,\ a_1=3,\ b_2=0$) was independently hand-verified during this audit by direct Lanczos iteration on the stated $3\times3$ matrix, and is correct — but the PNG itself is reproducible only by re-deriving a plotting script from scratch, which was not done here per the instruction not to invent code. See `SCIENTIFIC_ISSUES.md`.
- `blackstring_quiver/figures/` in the original project contained two additional PNGs (`02_J_independent_equilibrium.png`, `03_barrier_crossing_trajectories.png`) that the *current* version of `run_blackstring_quiver.py` does **not** regenerate (confirmed: their timestamps did not update on a clean re-run, while the three report-cited figures did). These are stale outputs of an earlier script version and were excluded from this repository; the originals are untouched in the source project. Not scientifically load-bearing (not cited by the report), but flagged for transparency.
- `Krylov_MIORPA.bib` (127 entries, in the original project root) is **not** referenced by any `.tex` file found (`\bibliography{Krylov}` in the final report uses `Krylov.bib`, 116 entries). Appears to be an earlier/duplicate bibliography from a draft. Left out of `report/` to avoid ambiguity about which bib file is authoritative; not deleted from the original project.
- The original project directory also contains presentation material (multiple `krylov_beamer*.tex` decks, `QA_prep_soutenance.tex`, `fiche_1page_reunion.tex`, `speech_notes_*.md`, `oxford_assets/`) and two earlier report drafts (`Report.tex`, `Plan.tex`) that are **not** included in this repository, since they are presentation/meeting artifacts rather than code or the final report. `Notes.tex` was the one exception included (as `report/supplementary/Notes.tex`) because it is directly cited by name in code comments as a derivation source.
- Two different Python virtual environments existed in the original project (`krylov_env`, Python 3.11.9, created on a different machine — its recorded interpreter path no longer resolves; `.venv`, Python 3.14.7). Neither is included in this repository (both are in `.gitignore`); `requirements.txt` instead pins the exact package versions verified during this audit.

## NOT VERIFIED

- `report/supplementary/Notes.tex` was not compiled or checked equation-by-equation against the final report during this audit; it is included only as a traceability reference for the derivations that code comments point to.
- The remaining 8 (of 11) figures in `quiver_comparison/figures/` beyond the 3 cited in the report (`02`–`08`, `09_f1_shape_vs_P`) were confirmed to regenerate byte-for-different-but-successfully (i.e. the scripts ran without error and produced a file at the expected path) but were not individually cross-checked against any report claim, since the report does not cite them.
- Cross-platform behavior (Linux/macOS) was not tested — all verification above was performed on Windows with TinyTeX and a Python 3.14 interpreter local to this machine.
- GitHub repository creation, first release/tag, and Zenodo archiving were **not** performed by this audit (cannot be done without the user's GitHub/Zenodo accounts) — see `ZENODO_RELEASE.md` for the exact manual steps remaining.

## ACTION REQUIRED (before or shortly after publishing)

1. Replace `REPLACE_WITH_GITHUB_USERNAME` in `README.md`, `CITATION.cff`, and `report/Holographic_Krylov_Complexity_MIORPA_Report.tex` with the actual GitHub username/organization once the repository is created.
2. Replace `REPLACE_WITH_ZENODO_DOI` / `TO_BE_REPLACED_AFTER_ZENODO_PUBLICATION` in the same three files after the Zenodo archive is minted (see `ZENODO_RELEASE.md`).
3. Add an ORCID to `CITATION.cff` if/when available — none was found in any project file, so none was invented.
4. Decide whether to keep `figures/ct_example.png` as a static asset with a documented "no reproducing script found" caveat (current state), or write a new small script to regenerate it — the latter was intentionally **not** done here since it would be new code, not recovered code.
5. Run `git init`, review `git status` for anything unexpected, then make the first commit and tag `v1.0.0` (exact commands in `README.md`/`CHANGELOG.md`).
