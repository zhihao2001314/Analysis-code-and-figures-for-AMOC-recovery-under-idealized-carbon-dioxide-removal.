# AMOC manuscript: local release review package

This is a draft for author review, not a published or licensed release. No GitHub upload or DOI registration has been performed.

## Contents

Four main figures and eight candidate supplementary figures are in `figures/`. Supplement numbering is provisional. `figure_index.csv` maps every figure to its plot, calculation and data directory. `output/` contains the selected processed data and their upstream dependencies. `src/` contains the analysis scripts and `plot/` the plotting scripts. Original input data are not bundled.

## Reproduce figures from processed results

Use Python and the packages in `requirements.txt`. From this directory run each distinct script in `plot_order.txt` in order. Plot scripts locate this package from their own file path. New renders go to `figure/`; the supplied reference images in `figures/` are retained. Before repeating a render, move its existing `figure/` directory to a dated archive rather than overwriting it. FigureS3–S5 are produced together by one plotting script.

## Recompute from original inputs

Analysis scripts are preserved with their scientific parameters and original local raw-data paths. Edit the path blocks for your machine before recomputing. Obtain the exact CMIP6 experiments, ensemble members, variables, grid versions and calendar windows in source inventories and metadata. Source-code dependencies and processed intermediates are retained, but a fully automated fresh ESGF download is not included. Do not run archive/revision scripts as raw-data calculations.

## Interpretation and known qualifications

Figure3 radiation uses two non-recovery models(CESM2,CNRM) and three recovery models; its temperature maps use three versus three including UKESM. MainFigure4 uses net precipitation(P−E), not precipitation alone. FigureS7 displays upward-positive evaporationE and satisfies Δ(P−E)=ΔP−ΔE. Group responses use unanimous three-model sign agreement; between-group differences use their own block-bootstrap/FDR fields. Figure1 andFigureS1 panelb is an idealized small-area closed schematic, not fitted observations.

The legacy regional/temperature pattern products retain their documented source windows, including GFDLinitial6–15 where specified. The new globaltas time series and newer atmospheric endpoint diagnostics use GFDL1–10. These are not identical initial-reference samples; this difference is explicitly preserved for review. CAS branch-time metadata and control-branch alignment have documented qualifications; these are not claimed to be resolved by packaging. See original metadata/readmes.

## Before publication

Confirm final supplementary numbering, article title, authors/affiliations, license and third-party redistribution permissions. Add the exact ESGF/IAP/NSIDC product citations and datasetDOIs, plus GitHub/Zenodo identifiers after release. No identifiers or licenses have been invented. Original global CMIP6 files remain at their source repositories; the derived arrays and tables used by the figures are supplied here. Large data should be distributed in the Zenodo data archive if unsuitable for regular GitHubgit.
