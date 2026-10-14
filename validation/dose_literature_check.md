# Dose-threshold prediction vs literature (full-text checks via Undermind, 2026-10-10)
- Les22 (diluted FMT, C. diff): graded loss of protection with dilution (colonized 6/14 undiluted -> 10/12 at 1:10 -> 10/14 at 1:10^2), antibiotic-background dependent; pathogen load fixed at 1e3 spores -> pathogen-load dependence UNTESTED.
- Hon23 (18-strain consortium): no dose titration; mechanism = gluconate competition (nutrient competition), cross-feeding not shown.
- Jen23 (B. longum cooperative vs B. producta competitive consortium, C. diff): single regimen, cooperative consortium survived; amino-acid cross-feeding (Bifidobacterium). No dose series.
Model (v4, 60 virtual hosts, lognormal CV 15% on 8 parameters, seed 20261010; stage2a/dose_heterogeneity_v4.py):
 P0=0.01: fraction protected 0 (dose<=0.06) -> 14/60 (0.1) -> 35/60 (0.15) -> 49/60 (0.2) -> 58/60 (0.3) -> 100% (>=0.5)
 P0=0.10: 0 (<=0.1) -> 2/60 (0.15) -> 9/60 (0.2) -> 30/60 (0.3) -> 53/60 (0.5) -> 100% (1.0)  => ED50 shifts ~0.14 -> ~0.3 with 10x pathogen load
Reading: graded fraction-vs-dose curve from host heterogeneity is qualitatively consistent with Les22; the pathogen-load shift is a novel testable prediction.
