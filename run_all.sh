#!/usr/bin/env bash
# Re-run every analysis and redraw the figures. Each script must be run from its own folder
# (they write CSV files next to themselves and read ../stage2a). Total time: roughly an hour on 2 cores.
set -e
cd "$(dirname "$0")"
( cd stage2a    && python analysis_stage2a_main_v4.py   | tee main_v4.log )
( cd stage2a    && python analysis_stage2a_robust_v4.py | tee robust_v4.log )
( cd stage2a    && python dose_heterogeneity_v4.py )
( cd stage2a    && python check_v4_gut.py )
( cd stage2b    && python design_matched_dose.py | tee design_matched_dose.log )
( cd validation && python rio15_representative.py | tee rio15_representative_output.txt )
( cd figures    && Rscript fig1_scheme.R && Rscript fig2_timecourse.R )
( cd stage2b    && Rscript fig_summary.R )
echo done
