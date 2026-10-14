# Testable prediction (matched total dose)
Design: fixed total inoculum (B + F), vary the B fraction 0..1, challenge with a fixed pathogen load.
Prediction (spatial and well-mixed models agree): single-species arms (fraction 0 or 1) never suppress the pathogen even at 4x the pair's threshold dose;
mixed arms do, over a broad ratio window (B fraction ~0.1-0.9 at high dose), narrowing to ~0.5-0.75 near the threshold. Total-dose threshold ~0.3-0.5 (well-mixed) / ~0.5 (spatial).
Falsification: a single species at matched total dose performs as well as the pair, or the pair's advantage disappears when acetate is supplied externally (control).
Spatial model: B concentrates in the lumen, F near the mucosa; pair threshold 0.18-0.30 per species over 5 diffusion settings (0.237 default); singles fail (inoculum up to 4); mucosa-seeded pathogen also cleared.
Python 3.13.16 numpy 2.5.3 scipy 1.18.1; R 4.3.3 ggplot2 3.4.4; deterministic (no seed).
