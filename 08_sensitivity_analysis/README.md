# Final Candidate Prioritization: Sensitivity Analysis

This directory contains reproducible post hoc sensitivity analyses
for nine quinone candidates targeting Madurella mycetomatis
scytalone dehydratase (MmSCD).

## Input
- ligand_metrics.csv: Six computational descriptors for nine ligands.

## Analyses
1. 01_equal_weight.py: Equal-weight composite scoring.
2. 02_leave_one_out.py: Leave-one-metric-out sensitivity analysis.
3. 03_rank_aggregation.py: Mean-rank aggregation.
4. 04_random_weights.py: Monte Carlo random-weight analysis.

## Figure scripts
- 05_heatmap.py
- 06_monte_carlo_plot.py
- 07_leave_one_out_bubble.py

## Reproduction

From the GitHub_Pipeline repository root, run:

python 08_sensitivity_analysis/01_equal_weight.py
python 08_sensitivity_analysis/02_leave_one_out.py
python 08_sensitivity_analysis/03_rank_aggregation.py
python 08_sensitivity_analysis/04_random_weights.py

Monte Carlo uses 100,000 iterations, six weights sampled
from Dirichlet(1,1,1,1,1,1), and random seed 42.

Expected result: 67,956/100,000 iterations (67.956%)
reproduce the complete original five-candidate shortlist.

## Software

Project environment inspected on Narval:
- NumPy 2.4.2+computecanada
- Pandas 3.0.5+computecanada
- Matplotlib 3.11.1+computecanada

Python 3.11.4 was observed in the Narval system environment.
The Python version within the project virtual environment
should be confirmed separately.

## Interpretation

These analyses evaluate the stability of computational
candidate prioritization, not experimental inhibitory activity.
