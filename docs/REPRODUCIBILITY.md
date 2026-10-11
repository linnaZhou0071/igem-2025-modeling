# Reproducing the models and analyses

## Requirements

The portable runner uses Python 3.10+, GNU Make and `g++` with OpenMP on Linux. It builds against the included PhysiCell 1.14.2 engine. A Studio GUI is not required. Original Makefiles use `-march=native`; compiled binaries should be rebuilt on each machine.

Scoring and offline HTML generation use the Python standard library. MAT analysis and regenerated static figures additionally need the packages in `requirements.txt`. Videos require a separate ffmpeg installation.

## Run one model

From the repository root:

```bash
python3 scripts/run_model.py Recurrence --threads 4
```

Valid model names are `AD_trial1`, `AD_trial2`, `Recurrence`, `Treatment`, `AD_3D` and `Recurrence_3D`. Results are written under `build/<model>/output/`. The first run compiles the engine; subsequent runs reuse that model's compiled objects.

```bash
# Generate working files without compiling or running.
python3 scripts/run_model.py AD_3D --prepare-only

# Compile only.
python3 scripts/run_model.py AD_3D --build-only --jobs 2

# A short run checks compilation, configuration loading and early dynamics.
python3 scripts/run_model.py AD_3D --threads 2 --max-time 1
```

The runner copies the original files to an isolated build folder, selects the active settings recorded in [RESEARCH.md](RESEARCH.md), and changes only input/output paths, OpenMP thread count, and optional duration. It uses each project's own `cells.csv` and `cell_rules.csv`. Biological parameters and archived population values remain intact. Published configuration and metadata copies have machine-specific paths removed, with source and exported hashes recorded in the manifest. For Recurrence, this resolves an original cross-reference to Treatment's rule directory; historical results may have used a different rule file.

Full historical configurations run for 2,880 minutes in most models, 5,760 in Treatment, and 2,881 in Recurrence_3D. A successful one-minute check verifies software execution, not the full research claims or long-term dynamics. Fresh stochastic trajectories and portability fixes need not reproduce a historical run bit-for-bit.

## Recompute archived scores and figures

```bash
python3 models/Recurrence/score_runs.py --output results/scored_runs.csv
python3 models/Recurrence/score_runs.py \
  --verify-run 791 --output /tmp/igem-scored-check.csv
python3 models/Recurrence/scatter3d_axes_offline.py \
  --csv results/scored_runs.csv --top-percent 0.05 \
  --output results/parameter-space.html

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/analyze_sweep.py
```

`analyze_sweep.py` checks timestamps, nonnegative finite counts, consistent parameter values, run IDs, and all six score components against the archive. It writes `results/sweep-verification.json` and three PNG plots under `docs/figures/`. It reuses the original scoring function, so the comparison checks archived-data consistency rather than independently validating the scoring implementation.

Each archived run includes 25 hourly population observations from 0 to 1,440 minutes. This is a **24-hour archived analysis**, even though the supplied active recurrence configuration specifies 48 hours. One trajectory is supplied per grid point, without replicate-based uncertainty estimates.

## Start a new sweep

The original `parameter_sweep.py` expects a traditional layout with `user_projects/Recurrence` and `bin/Recurrence`. It is retained for inspection, but should not be launched in place from `models/Recurrence`: its paths reflect the original installation and it modifies the working rule file.

To use it in an isolated working installation:

```bash
python3 scripts/run_model.py Recurrence --build-only
python3 - <<'PY'
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
work = Path('build/Recurrence').resolve()
project = work / 'user_projects/Recurrence'
shutil.copytree('models/Recurrence', project, dirs_exist_ok=True,
                ignore=shutil.ignore_patterns('sweep_outputs', '__pycache__'))
shutil.copy2(work / 'config/portable.xml', project / 'config/settings.xml')
tree = ET.parse(project / 'config/settings.xml')
settings = tree.getroot()
for parent in (settings.find('initial_conditions/cell_positions'),
               *settings.findall('cell_rules/rulesets/ruleset')):
    if parent is not None:
        parent.find('folder').text = str(project / 'config')
settings.find('save/folder').text = str(work / 'output')
tree.write(project / 'config/settings.xml', encoding='utf-8', xml_declaration=True)
(work / 'bin').mkdir(exist_ok=True)
shutil.copy2(work / 'project', work / 'bin/Recurrence')
PY
cd build/Recurrence
# Run a single parameter combination before considering all 1,000.
python user_projects/Recurrence/parameter_sweep.py --test
```

The historical sweep restores `cell_rules_backup.csv` before modifying two secretion saturations. That backup is omitted from the curated source; the script creates a backup from the supplied active rules in its new working copy. New runs therefore use the curated active rules, which are not guaranteed to match the original historical sweep baseline. Archived `metadata.json` files retain historical timestamps and parameters; machine paths have been replaced with relative CSV filenames for publication.

## Inspect fresh spatial output

After a sufficiently long simulation, inspect mesh dimensions and substrate names in the generated XML before choosing `--grid-shape` and substrate indices:

```bash
python models/AD_trial2/tools/extract_population_timeseries.py \
  --output-dir build/AD_trial2/output \
  --out-csv build/AD_trial2/population_timeseries.csv

python models/AD_trial2/tools/visualize_physicell_mat.py \
  build/AD_trial2/output --substrates 3 7 \
  --grid-shape 50 50 1 --outdir build/AD_trial2/viz --plot
```

The numeric indices and grid shape above are examples; confirm them for the chosen model. The original visualization tools also support substrate exports, videos and bivariate Moran's I; see their CLI help and included tool documentation.

For a separately installed PhysiCell Studio, `-c` selects the configuration and `-e` the executable. Studio is optional and is not bundled in this repository.

## Availability limits

- No executable MATLAB/SimBiology source was included, so the ODE chapter cannot be rerun from these files alone.
- Full historical MAT/XML spatial trajectories and all referenced simulation movies were not supplied; embedded wiki figures/animations are preserved.
- The supplied Recurrence_3D rule table and custom reaction function do not contain a direct VHH-to-cytokine suppression term. The supplied implementation therefore does not, by itself, substantiate every VHH-related claim in the historical 3D results draft. No new biological rule has been invented to fill this gap.
- Some original 3D generator scripts use hard-coded paths and mutate CSVs. Inspect them and work on a copy; the portable runner starts from the supplied final CSVs.

Verification results and scope are recorded in `results/model-smoke-checks.json` and `results/sweep-verification.json`.
