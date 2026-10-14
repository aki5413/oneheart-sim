# oneheart-sim

Code, output data and figures for the manuscript

> Yamada A. *A cross-feeding pair of gut bacteria creates a dose threshold for pathogen exclusion: a mechanistic model.* (submitted)

The model is a chemostat consumer-resource model in which *Bifidobacterium* (B) releases acetate and breakdown products,
*Faecalibacterium* (F) turns acetate into butyrate, and butyrate-fuelled host oxygen use removes the oxygen a pathogen (P) needs.
Concentrations and biomass are in **arbitrary units** unless a run is compared with published millimolar data.

## What is where

| Folder | Content |
|---|---|
| `stage2a/` | Well-mixed model. `crossfeed_v4.py` is the model (state: sugar C, host carbon H, acetate, butyrate, oxygen O, B, F, P, breakdown pool M). `crossfeed_v3.py` is the earlier version (v4 with `lk=0` or `wM=0` equals v3). Analysis scripts `*_v4.py` and their CSV outputs `*_v4.csv`. Files without `_v4` are the earlier v3 results, kept for the record. |
| `stage2b/` | One-dimensional lumen-to-mucosa model (`spatial_v4.py`) and the matched-total-dose design (`design_matched_dose.py`). |
| `validation/` | Check against the B + F co-culture of Ríos-Covián et al. 2015 (FEMS Microbiol Lett, doi:10.1093/femsle/fnv176). |
| `figures/`, `stage2b/fig_summary.R` | R scripts and PDF/PNG for Fig. 1-3 of the manuscript. |

## Which script gives which result

| Manuscript result | Script | Output |
|---|---|---|
| Two stable states, time courses (Fig. 2) | `stage2a/analysis_stage2a_main_v4.py` | `s2a_timecourse_v4.csv`, `s2a_equilibria_v4.csv` |
| Pair threshold 0.146; critical cross-feeding 3.59 and 9.97 (Fig. 3a) | same | `s2a_threshold_vs_aB_v4.csv`, `s2a_critical_v4.csv`, `s2a_invasion_v4.csv` |
| Ablations, validity map | same | `s2a_ablation_v4.csv`, `s2a_boundary_map_v4.csv` |
| Threshold vs pathogen load (Fig. 3b) | same | `s2a_threshold_vs_P0_v4.csv` |
| 300 random parameter sets, Spearman correlations | `stage2a/analysis_stage2a_robust_v4.py` | `s2a_robustness_v4.csv`, `s2a_sensitivity_v4.csv` |
| Host variation, graded response (Fig. 3c) | `stage2a/dose_heterogeneity_v4.py` | `s2a_dose_heterogeneity_v4.csv` |
| Gut-mode threshold with and without the breakdown-product pathway | `stage2a/check_v4_gut.py` | printed output |
| Spatial threshold 0.237; diffusion sensitivity 0.18-0.30; matched-dose design (Fig. 3d) | `stage2b/design_matched_dose.py` | `s2b_matched_dose.csv`, `design_matched_dose.log` |
| Co-culture check (butyrate, acetate) | `validation/rio15_representative.py` | `rio15_representative_output.txt` |
| Parameter scan that led to the v4 breakdown-product pathway | `validation/scan_rio15_v4.py`, `scan_rio15_v4b.py` | `rio15_scan_v4*.csv`, `validation_summary_v4.txt` |

## How to run

```bash
pip install -r requirements.txt          # Python 3.13; numpy, scipy, pandas
bash run_all.sh                          # everything (about an hour on 2 cores)
```

Or run one script from inside its own folder, because scripts write their CSV next to themselves and import from `../stage2a`:

```bash
cd stage2a && python analysis_stage2a_main_v4.py
cd ../validation && python rio15_representative.py     # seconds
```

R packages for the figures are listed in `R_packages.txt`.

## Reproducibility

- Differential equations are integrated with LSODA (SciPy), relative tolerance 1e-7 to 1e-8.
- Random numbers: `numpy.random.default_rng(20261010)` in the robustness and host-variation scripts. Threshold searches (bisection, 14 steps, range 0.001-10) and the spatial runs are deterministic.
- Versions used: Python 3.13.16 (numpy 2.5.3, scipy 1.18.1, pandas 3.0.5) for the v4 analyses; Python 3.12.3 (numpy 2.4.4, scipy 1.17.1, pandas 3.0.2) for the earlier v3 runs. R 4.3.3, ggplot2 3.4.4, patchwork 1.2.0. Each analysis writes an `environment_*.txt` file.
- Small numerical differences (last digits) can appear with other library versions.

## Limitations to keep in mind

- Units are arbitrary; doses cannot be converted to colony-forming units.
- The co-culture reference values come from the text of Ríos-Covián et al. 2015, not from digitised figures. The breakdown-product pathway was tuned to them.
- The spatial model is one-dimensional with gradients chosen by the author.
- This is a hypothesis-generating model and has not been tested in vivo.

## Use of AI tools

The code, figures and first drafts of the text were produced with the large language model Claude (Anthropic) under the author's direction.
The author designed the study and is responsible for the content.

## Licence and citation

Code: MIT (`LICENSE`). Data files (CSV) and figures: CC BY 4.0. Please cite the manuscript and the Zenodo record (`CITATION.cff`).
