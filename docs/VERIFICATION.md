# Verification appendix — October 2026

This appendix documents work added during portfolio preparation. The following plots were newly generated from the archived 2025 CSV files in October 2026. They are not original 2025 wiki figures and do not represent new simulation experiments.

The [original 2025 results figures](wiki/physicell4.md) remain the visual record of the iGEM project.

## Archived-data consistency

All 1,000 archived runs were rescored using the original scoring function. Each contains 25 observations from 0 to 1,440 minutes. This is the 24-hour archived sweep, distinct from the current standalone `models/Recurrence/config/settings.xml` configuration (`max_time = 2,880` minutes, or 48 hours). Rescoring the archive does not rerun a 48-hour simulation. All six penalty components and final scores agree exactly with the archived score table. This is a data consistency check, not independent biological validation.

![2026 redraw of archived parameter scores](figures/parameter-space.png)

*2026 redraw: orange rings mark the 50 highest-scoring grid points.*

![2026 comparison of highest- and lowest-scoring archived trajectories](figures/population-comparison.png)

*2026 comparison: the highest- and lowest-scoring archived trajectories, selected by the original population objective. These particular comparative panels were not present in the original wiki.*

![2026 mean-score summary across AIP-C settings](figures/mean-score.png)

*2026 summary: mean score across the ten AIP-C settings for each VHH parameter pair.*

## Independent heatmap verification — 15 October 2026

The heatmap above was checked against all 1,000 archived population time series using a separate implementation of the documented equations, without importing the original scorer or plotting script. The largest difference across scores and all six penalty components was **2.22 × 10⁻¹⁶**. Both archived score tables agree with that calculation.

The archive contains a complete 10 × 10 × 10 grid. Each of the 100 displayed cells is the arithmetic mean of **exactly ten scores**, one per AIP-C secretion setting, with the VHH–cytokine coefficient and VHH secretion saturation held fixed:

```math
\overline{S}(d,v) = \frac{1}{10}\sum_{j=1}^{10} S(d,a_j,v).
```

The x-axis increases from 0.1 to 1.0; the y-axis increases from 0.01 at the bottom to 0.5 at the top. All 100 cell-center PNG colors match their independently calculated means exactly (maximum RGB channel difference: **0/255**). The existing image was retained.

The largest cell mean is **0.906959**, at coefficient **0.8** and VHH secretion saturation **0.0644444**. This is an average over ten AIP-C settings; the highest individual run has score **0.994872**. Averaging across parameter settings does not provide a replicate-based uncertainty estimate or validate biological efficacy.

Inspect the [100-cell numerical table](../results/mean-score-grid.csv), [machine-readable verification](../results/mean-score-verification.json), and [independent checking script](../scripts/verify_sweep_summary.py). After installing `requirements.txt`, reproduce these checks with:

```bash
python scripts/verify_sweep_summary.py --check-image
```

## Execution checks

Six PhysiCell models compiled and completed one-minute runs using portable settings. This establishes configuration loading and early software execution, not full-horizon reproduction of the research findings.

Evidence: [model smoke checks](../results/model-smoke-checks.json) and [sweep consistency report](../results/sweep-verification.json). Regenerate this appendix’s plots with `python scripts/analyze_sweep.py` after installing `requirements.txt`.
