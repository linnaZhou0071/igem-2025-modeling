
import xml.etree.ElementTree as ET
import numpy as np
import subprocess
import os
try:
    import pandas as pd  # optional; we can fall back to plain CSV writes
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

try:
    import scipy.io as sio
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
import shutil
from pathlib import Path

# --- Configuration ---
# Resolve important paths (absolute to avoid CWD issues)
PROJECT_DIR = Path(__file__).resolve().parent
ROOT_DIR = PROJECT_DIR.parents[1]  # PhysiCell root
BIN_DIR = ROOT_DIR / 'bin'

# Base PhysiCell settings file (original provided by user)
CONFIG_FILE = PROJECT_DIR / 'config' / 'settings.xml'
"""We will NOT create or rely on a persistent modified working copy anymore.
Per run we will create a temporary copy that only changes the single user parameter
`vHH_to_cytokine_degradation_rate` (and nothing else in settings.xml) in accordance
with the user's instruction to modify only the specified parameter location.
"""

# Cell rules file paths
CELL_RULES_FILE = PROJECT_DIR / 'config' / 'cell_rules.csv'
CELL_RULES_BACKUP = PROJECT_DIR / 'config' / 'cell_rules_backup.csv'

# PhysiCell executable
EXECUTABLE = BIN_DIR / 'Recurrence'

# Directory to store all run outputs
BASE_OUTPUT_DIR = PROJECT_DIR / 'sweep_outputs'
# File to store summarized population data
SUMMARY_FILE = PROJECT_DIR / 'population_data.csv'

# --- Parameter Ranges ---
# vHH_to_cytokine_degradation_rate: (0, 1]
p1_degradation_rates = np.linspace(0.1, 1.0, 10)

# AIP_sa increases AIP_C secretion saturation: [0.5, 3]
p2_aipc_saturation_values = np.linspace(0.5, 3.0, 10)

# AIP_sa increases vHH secretion saturation: [0.01, 0.5]
p3_vhh_saturation_values = np.linspace(0.01, 0.5, 10)

# NOTE: We no longer override save interval, threads, or max_time inside settings.xml.
# The sweep strictly modifies ONLY:
#   1. user_parameters/vHH_to_cytokine_degradation_rate
#   2. In cell_rules.csv: the saturation (5th column) for the two lines:
#        engr_S.epidermis,AIP_sa,increases,AIP_C secretion
#        engr_S.epidermis,AIP_sa,increases,vHH secretion
# All other configuration values remain exactly as in the provided files.

# Retain constants for logic (e.g. deriving fallback times if needed)
SAVE_INTERVAL_MIN = 60  # used only for fallback time derivation if XML time not parsed

# Optional pyMCDS import (Studio logic). We'll only use when user passes --use-pymcds.
HAS_PYMCDS = False
PYMCDS_IMPORT_ERR = None
PYMCDS_SUPPORTS_MICROENV = True
STUDIO_BIN = ROOT_DIR / 'studio' / 'bin'
if STUDIO_BIN.exists():
    import sys
    if str(STUDIO_BIN) not in sys.path:
        sys.path.insert(0, str(STUDIO_BIN))
    try:
        from pyMCDS import pyMCDS  # type: ignore
        HAS_PYMCDS = True
        # Detect whether signature includes 'microenv'
        try:
            import inspect as _inspect
            sig = _inspect.signature(pyMCDS.__init__)
            if 'microenv' not in sig.parameters:
                PYMCDS_SUPPORTS_MICROENV = False
        except Exception:  # noqa: BLE001
            pass
    except Exception as e:  # noqa: BLE001
        PYMCDS_IMPORT_ERR = str(e)

def _read_base_output_folder():
    """Parse the base settings.xml to discover the (unmodified) output folder path.
    We need this so we can move generated output files into run-specific folders
    after each simulation without altering the settings file's folder node.
    Returns absolute Path or None if not found.
    """
    try:
        tree = ET.parse(CONFIG_FILE)
        root = tree.getroot()
        folder_node = root.find('.//save/folder')
        if folder_node is None or not folder_node.text:
            return None
        # Resolve relative path w.r.t. project root (same directory as CONFIG_FILE parent)
        base_dir = CONFIG_FILE.parent.parent  # project_dir/config -> project_dir
        rel_path = folder_node.text.strip()
        candidate = (base_dir / rel_path)
        if not candidate.exists():
            # Try resolving from repository root instead (PhysiCell root)
            repo_root = ROOT_DIR
            alt = (repo_root / rel_path)
            if alt.exists():
                return alt.resolve()
        return candidate.resolve()
    except Exception as e:  # noqa: BLE001
        print(f"[warn] Could not read base output folder from settings.xml: {e}")
        return None

def modify_only_parameters_and_rules(run_id, p1, p2, p3, disable_svg=True):
    """Create a temp copy of settings.xml, modifying ONLY vHH_to_cytokine_degradation_rate and (optionally)
    disabling SVG output (by setting <save><SVG><enable>false</enable></SVG> in the *temporary* file only).
    Also adjust only the two engr_S.epidermis rule lines in cell_rules.csv (5th column).
    Return (temp_config_path, run_specific_output_folder_path).
    """
    # Prepare per-run destination (where we'll MOVE outputs afterward)
    run_output_folder = BASE_OUTPUT_DIR / f'run_{run_id:03d}'
    run_output_folder.mkdir(parents=True, exist_ok=True)

    # Create temp config with ONLY the user parameter altered
    tree = ET.parse(CONFIG_FILE)
    root = tree.getroot()
    p1_node = root.find('.//user_parameters/vHH_to_cytokine_degradation_rate')
    if p1_node is not None:
        p1_node.text = f'{p1:.6g}'
    else:
        print("Warning: vHH_to_cytokine_degradation_rate not found in user_parameters")

    # Optionally disable SVG output ONLY in this temporary config (does not touch original file)
    if disable_svg:
        svg_enable_node = root.find('.//save/SVG/enable')
        if svg_enable_node is not None:
            orig_val = svg_enable_node.text.strip().lower() if svg_enable_node.text else 'true'
            if orig_val != 'false':
                svg_enable_node.text = 'false'
        else:
            print('[warn] Could not find <save><SVG><enable> node to disable SVG output.')
    temp_config_path = PROJECT_DIR / 'config' / f'temp_settings_run_{run_id}.xml'
    tree.write(temp_config_path)

    # Adjust cell_rules.csv starting from pristine backup
    if not CELL_RULES_BACKUP.exists():
        shutil.copy2(CELL_RULES_FILE, CELL_RULES_BACKUP)
    shutil.copy2(CELL_RULES_BACKUP, CELL_RULES_FILE)

    lines = CELL_RULES_FILE.read_text().strip().splitlines()
    new_lines = []
    for line in lines:
        parts = line.split(',')
        if len(parts) >= 8:
            cell, signal, direction, behavior = parts[0], parts[1], parts[2], parts[3]
            if cell == 'engr_S.epidermis' and signal == 'AIP_sa' and direction == 'increases':
                if behavior == 'AIP_C secretion':
                    parts[4] = f'{p2:.6g}'
                elif behavior == 'vHH secretion':
                    parts[4] = f'{p3:.6g}'
                line = ','.join(parts)
        new_lines.append(line)
    CELL_RULES_FILE.write_text('\n'.join(new_lines) + '\n')

    return str(temp_config_path), str(run_output_folder)

def run_simulation(run_id, config_path):
    """
    Runs the PhysiCell simulation for a given configuration.
    """
    print(f"--- Starting Run {run_id:03d} ---")
    command = [str(EXECUTABLE), config_path]
    try:
        # Using subprocess.run to wait for completion
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout.splitlines()[0][:200])
        print(f"--- Finished Run {run_id:03d} ---")
        return True
    except subprocess.CalledProcessError as e:
        print(f"!!! Error in Run {run_id:03d} !!!")
        print(f"Command: {' '.join(command)}")
        print(f"Return code: {e.returncode}")
        print(f"Output:\n{e.stdout}")
        print(f"Error Output:\n{e.stderr}")
        return False
    finally:
        # Keep temp config until after we relocate outputs (handled outside) then delete.
        pass

def relocate_outputs_to_run_folder(base_output_folder: Path, run_output_folder: Path):
    """Move newly generated output files from the (constant) base output folder to the
    per-run folder. We include typical patterns. If folders differ, perform moves.
    """
    if base_output_folder.resolve() == run_output_folder.resolve():
        return  # nothing to do (shouldn't happen in our design)
    import glob
    patterns = [
        'output*.xml',
        'output*_cells.mat',
        'output*.csv',
        '*.svg',  # SVGs should no longer be generated if disable_svg=True, but keep pattern for safety
        'initial.svg'
    ]
    for pat in patterns:
        for src in glob.glob(str(base_output_folder / pat)):
            dst = run_output_folder / os.path.basename(src)
            try:
                shutil.move(src, dst)
            except Exception as e:  # noqa: BLE001
                print(f"[warn] Failed to move {src} -> {dst}: {e}")

def clean_base_output_folder(base_output_folder: Path):
    """Remove prior outputs so each simulation produces a clean set before relocation."""
    import glob
    patterns = [
        'output*.xml', 'output*_cells.mat', 'output*.csv', '*.svg', 'initial.svg'
    ]
    for pat in patterns:
        for f in glob.glob(str(base_output_folder / pat)):
            try:
                os.remove(f)
            except OSError:
                pass

def purge_run_directories(keep_files=('population_timeseries.csv',), verbose=True):
    """Delete every file inside each run_* directory under BASE_OUTPUT_DIR except the listed keep_files.
    Returns (num_dirs_processed, num_files_deleted).
    """
    if not BASE_OUTPUT_DIR.exists():
        if verbose:
            print(f"[purge] Base sweep output dir {BASE_OUTPUT_DIR} does not exist; nothing to purge.")
        return 0,0
    dirs = sorted([d for d in BASE_OUTPUT_DIR.iterdir() if d.is_dir() and d.name.startswith('run_')])
    deleted = 0
    for d in dirs:
        for item in d.iterdir():
            if item.is_file() and item.name not in keep_files:
                try:
                    item.unlink()
                    deleted += 1
                except OSError as e:
                    if verbose:
                        print(f"[purge] Failed to delete {item}: {e}")
            elif item.is_dir():
                # remove whole subdirectory tree
                try:
                    shutil.rmtree(item)
                    deleted += 1
                except OSError as e:
                    if verbose:
                        print(f"[purge] Failed to delete directory {item}: {e}")
    if verbose:
        print(f"[purge] Processed {len(dirs)} run directories. Deleted {deleted} entries.")
    return len(dirs), deleted

def extract_population_data(output_folder, run_id, p1, p2, p3, include_dead=False):
    """Extract population time series by reading outputXXXX.xml for time and corresponding *_cells.mat for counts.

            Logic (MAT fallback): For each output%08d.xml, load output%08d_cells.mat.
            Labels: index 5 = cell_type, index 6 = cycle_model, index 26 = dead.
            To match Studio population plot we should use cycle_model < 100 unless include_dead requested.
    """
    if not HAS_SCIPY:
        print("Warning: scipy not installed; cannot extract populations.")
        return None

    xml_files = sorted([f for f in os.listdir(output_folder) if f.startswith('output') and f.endswith('.xml')])
    if not xml_files:
        print(f"Warning: no output*.xml files in {output_folder}")
        return None

    times = []
    s_counts = []
    engr_counts = []

    import re, numpy as np
    time_pattern = re.compile(r'<current_time[^>]*>([-+0-9.eE]+)')

    for xf in xml_files:
        base = xf.replace('.xml','')
        mat_file = os.path.join(output_folder, base + '_cells.mat')
        xml_path = os.path.join(output_folder, xf)
        if not os.path.exists(mat_file):
            continue
        # parse time
        t_val = None
        with open(xml_path,'r') as f:
            for line in f:
                m = time_pattern.search(line)
                if m:
                    try:
                        t_val = float(m.group(1))
                    except:
                        t_val = None
                    break
        if t_val is None:
            # derive from index * SAVE_INTERVAL_MIN
            try:
                idx = int(base.replace('output',''))
                t_val = idx * SAVE_INTERVAL_MIN
            except:
                t_val = 0.0
        # load .mat
        try:
            mat = sio.loadmat(mat_file)
        except Exception as e:
            print(f"  MAT load failure {mat_file}: {e}")
            continue
        # MAT layout here is (features, n_cells); earlier assumption was transposed.
        cells_matrix = None
        for k,v in mat.items():
            if k.startswith('__'):
                continue
            if isinstance(v, np.ndarray) and v.ndim==2 and v.shape[0] < v.shape[1]:
                # heuristic: features (rows) less than cells (cols)
                cells_matrix = v
                break
        if cells_matrix is None:
            continue
        # Identify cell_type row (contains integers 0..5) and cycle row (few distinct states  (e.g., 4)).
        # Cache detection once.
        if 'EXTRACT_CACHE' not in locals():
            EXTRACT_CACHE = {}
        if 'cell_type_row' not in EXTRACT_CACHE:
            import math
            ct_row = None; cyc_row = None
            for ri in range(cells_matrix.shape[0]):
                u = np.unique(cells_matrix[ri,:])
                if 2 <= len(u) <= 10 and np.allclose(u, np.round(u)) and set(u).issuperset({0,5}):
                    ct_row = ri
                if cyc_row is None and 2 <= len(u) <= 10 and np.allclose(u, np.round(u)) and not set(u).issuperset({0,5}):
                    cyc_row = ri
            # Fallback defaults from current data inspection
            if ct_row is None:
                ct_row = 5
            if cyc_row is None:
                cyc_row = 4
            EXTRACT_CACHE['cell_type_row'] = ct_row
            EXTRACT_CACHE['cycle_row'] = cyc_row
        ct_row = EXTRACT_CACHE['cell_type_row']
        cyc_row = EXTRACT_CACHE['cycle_row']
        cell_types = cells_matrix[ct_row,:].astype(int)
        cycle_model = cells_matrix[cyc_row,:]
        if include_dead:
            live_mask = np.ones_like(cycle_model, dtype=bool)
        else:
            # For now we still use cycle_model < 100; if cycle_model row is phases (0..3) all are live.
            live_mask = (cycle_model < 100.)
        s_ct = int(np.count_nonzero( (cell_types == 0) & live_mask ))
        e_ct = int(np.count_nonzero( (cell_types == 5) & live_mask ))
        times.append(t_val)
        s_counts.append(s_ct)
        engr_counts.append(e_ct)

    if not times:
        return None

    order = sorted(range(len(times)), key=lambda i: times[i])
    times = [times[i] for i in order]
    s_counts = [s_counts[i] for i in order]
    engr_counts = [engr_counts[i] for i in order]

    # Always write per-run CSV
    per_run_csv = BASE_OUTPUT_DIR / f'run_{run_id:03d}' / 'population_timeseries.csv'
    with open(per_run_csv,'w') as f:
        f.write('time_min,S_aurenus_count,engr_S_epidermis_count,run_id,degradation_rate,aipc_saturation,vhh_saturation\n')
        for t,s,e in zip(times,s_counts,engr_counts):
            f.write(f'{t},{s},{e},{run_id},{p1},{p2},{p3}\n')

    # Always return a simple dict so caller can append even without pandas.
    result = { 'run_id': run_id, 'p1': p1, 'p2': p2, 'p3': p3, 'csv': str(per_run_csv) }
    return result

def extract_population_data_pymcds(output_folder, run_id, p1, p2, p3, include_dead=False, microenv=False, exact_times=False):
    """Replicate Studio (vis_base.py cell_counts_cb) logic using pyMCDS.

    Steps:
      1. Collect sorted output*.xml basenames.
      2. For each: mcds = pyMCDS(basename, output_dir, microenv=False, graph=False, verbose=False)
      3. time vector: np.linspace(0, last_time, len(files))  (exactly what Studio does)
      4. counts: for each snapshot, number of cells with cell_type==0 (S.aurenus) & cycle_model<100
                 and cell_type==5 (engr_S.epidermis) & cycle_model<100 (unless include_dead)
    Writes per-run CSV identical header to MAT method.
    """
    import glob
    if not HAS_PYMCDS:
        print(f"pyMCDS unavailable: {PYMCDS_IMPORT_ERR}. Falling back to MAT method.")
        return extract_population_data(output_folder, run_id, p1, p2, p3, include_dead=include_dead)

    pattern = os.path.join(output_folder, 'output*.xml')
    xml_files = sorted(glob.glob(pattern))
    if not xml_files:
        print(f"Warning: no output*.xml files in {output_folder}")
        return None

    snapshots = []
    for xf in xml_files:
        basename = os.path.basename(xf)
        try:
            if PYMCDS_SUPPORTS_MICROENV:
                if microenv:
                    m = pyMCDS(basename, output_folder, microenv=True, graph=False, verbose=False)
                else:
                    m = pyMCDS(basename, output_folder, microenv=False, graph=False, verbose=False)
            else:
                # Older pyMCDS version without microenv arg
                if microenv:
                    print("[pyMCDS] Installed version lacks 'microenv' parameter; proceeding without microenvironment data.")
                m = pyMCDS(basename, output_folder)  # type: ignore
            snapshots.append(m)
        except Exception as e:  # noqa: BLE001
            print(f" pyMCDS load failed for {basename}: {e}")
            continue
    if not snapshots:
        return None

    import numpy as _np
    if exact_times:
        times = _np.array([m.get_time() for m in snapshots], dtype=float)
    else:
        last_time = snapshots[-1].get_time()
        times = _np.linspace(0, last_time, len(snapshots))

    s_counts = []
    engr_counts = []
    for m in snapshots:
        data = m.data['discrete_cells']['data']
        ct = data['cell_type']
        # Some versions may store cycle model phases in a different key or not at all.
        cyc = data['cycle_model'] if 'cycle_model' in data else _np.zeros_like(ct)
        if include_dead:
            live_mask = _np.ones_like(ct, dtype=bool)
        else:
            live_mask = (cyc < 100.)
        s_counts.append(int(_np.count_nonzero( (ct == 0) & live_mask )))
        engr_counts.append(int(_np.count_nonzero( (ct == 5) & live_mask )))

    per_run_csv = BASE_OUTPUT_DIR / f'run_{run_id:03d}' / 'population_timeseries.csv'
    with open(per_run_csv,'w') as f:
        f.write('time_min,S_aurenus_count,engr_S_epidermis_count,run_id,degradation_rate,aipc_saturation,vhh_saturation\n')
        for t,s,e in zip(times,s_counts,engr_counts):
            f.write(f'{t},{s},{e},{run_id},{p1},{p2},{p3}\n')

    result = { 'run_id': run_id, 'p1': p1, 'p2': p2, 'p3': p3, 'csv': str(per_run_csv) }
    return result

def main(test_mode=False, aggregate=False, include_dead=False, use_pymcds=False, microenv=False, exact_times=False,
         resume=True, redo_failed=False, start_run=None, end_run=None, keep_only_population=False, disable_svg=True):
    """Main orchestration with resume/checkpoint support."""
    print("Starting parameter sweep...")
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    param_combinations = [(p1, p2, p3) for p1 in p1_degradation_rates for p2 in p2_aipc_saturation_values for p3 in p3_vhh_saturation_values]
    if test_mode:
        print("\n*** RUNNING IN TEST MODE ***")
        print("This will only perform ONE simulation run.\n")
        param_combinations = param_combinations[:1]

    total_runs = len(param_combinations)
    all_population_data = []  # when pandas available these become DataFrames; else dicts

    # Resume bookkeeping
    completed = set()
    param_grid_file = BASE_OUTPUT_DIR / 'param_grid.csv'
    if resume and BASE_OUTPUT_DIR.exists():
        for d in sorted(BASE_OUTPUT_DIR.glob('run_*')):
            rid_str = d.name.replace('run_','')
            if rid_str.isdigit():
                rid = int(rid_str)
                if (d / 'population_timeseries.csv').exists():
                    completed.add(rid)
        if completed:
            print(f"[resume] Found {len(completed)} completed runs.")

    # Write / reuse param grid
    if not param_grid_file.exists():
        with open(param_grid_file,'w') as fpg:
            fpg.write('run_id,degradation_rate,aipc_saturation,vhh_saturation\n')
            for rid,(p1,p2,p3) in enumerate(param_combinations):
                fpg.write(f"{rid},{p1},{p2},{p3}\n")
        print(f"[grid] Wrote {param_grid_file}")
    else:
        print(f"[grid] Using existing {param_grid_file}")

    runs_skipped = 0
    runs_executed = 0

    # Discover constant base output folder from config once
    base_output_folder = _read_base_output_folder()
    if base_output_folder is None:
        print('[error] Could not determine base output folder from settings.xml; aborting.')
        return

    for idx,(p1,p2,p3) in enumerate(param_combinations):
        run_id = idx
        if start_run is not None and run_id < start_run:
            continue
        if end_run is not None and run_id > end_run:
            break

        run_dir = BASE_OUTPUT_DIR / f'run_{run_id:03d}'
        pop_file = run_dir / 'population_timeseries.csv'

        if resume and (run_id in completed) and not redo_failed:
            print(f"[skip] run_{run_id:03d} already complete.")
            runs_skipped += 1
            continue
        if run_dir.exists() and pop_file.exists() and not redo_failed:
            print(f"[skip] run_{run_id:03d} population_timeseries.csv present.")
            runs_skipped += 1
            continue

        print(f"\nProcessing combination {idx+1}/{total_runs} (run_id={run_id:03d})")
        print(f"Parameters: degradation_rate={p1:.4f}, aipc_sat={p2:.4f}, vhh_sat={p3:.4f}")

        temp_config, run_specific_folder = modify_only_parameters_and_rules(run_id, p1, p2, p3, disable_svg=disable_svg)

        # Ensure base output folder clean before run
        clean_base_output_folder(base_output_folder)

        success = run_simulation(run_id, temp_config)

        # Move outputs produced in base folder into run-specific folder
        relocate_outputs_to_run_folder(base_output_folder, Path(run_specific_folder))

        # Remove the temp config file
        try:
            os.remove(temp_config)
        except OSError:
            pass

        if success:
            if use_pymcds:
                population_data = extract_population_data_pymcds(run_specific_folder, run_id, p1, p2, p3, include_dead=include_dead, microenv=microenv, exact_times=exact_times)
            else:
                population_data = extract_population_data(run_specific_folder, run_id, p1, p2, p3, include_dead=include_dead)
            if population_data is not None:
                all_population_data.append(population_data)
                # Write metadata json (parameters + death definition) in run folder
                try:
                    import json, time as _time
                    meta = {
                        'run_id': run_id,
                        'parameters': {
                            'vHH_to_cytokine_degradation_rate': p1,
                            'AIP_C_saturation': p2,
                            'vHH_saturation': p3
                        },
                        'death_definition': 'cycle_model >= 100 considered dead (Studio parity)',
                        'include_dead_counts': bool(include_dead),
                        'population_timeseries': population_data.get('csv') if isinstance(population_data, dict) else 'population_timeseries.csv',
                        'timestamp': _time.strftime('%Y-%m-%d %H:%M:%S')
                    }
                    with open(Path(run_specific_folder)/'metadata.json','w') as jf:
                        json.dump(meta, jf, indent=2)
                except Exception as e:  # noqa: BLE001
                    print(f"[warn] could not write metadata.json: {e}")
                if keep_only_population:
                    # Purge everything else in this run directory now
                    purge_run_directories(verbose=False)
        else:
            print(f"[warn] Run {run_id:03d} failed; no population file.")

        if aggregate and HAS_PANDAS and (run_id + 1) % 10 == 0 and all_population_data:
            print(f"[aggregate] Saving intermediate results to {SUMMARY_FILE}...")
            pd.concat(all_population_data, ignore_index=True).to_csv(SUMMARY_FILE, index=False)

        runs_executed += 1

    if aggregate and HAS_PANDAS and all_population_data:
        # Filter to DataFrames only
        dfs = [d for d in all_population_data if HAS_PANDAS and hasattr(d, 'to_csv')]
        if dfs:
            print(f"\n[aggregate] Saving all data to {SUMMARY_FILE}...")
            pd.concat(dfs, ignore_index=True).to_csv(SUMMARY_FILE, index=False)
            print("Done.")
    print("All per-run time series saved under each run_xxx folder.")

    print(f"\n[summary] Executed: {runs_executed}, Skipped: {runs_skipped}, Total combos: {total_runs}")

def analyze_initial_state():
    """Diagnose discrepancy between expected thousands of initial cells and observed few hundred.
    Reports:
      1. cells.csv total non-header lines and a sample of first/last 3 type labels.
      2. Distinct type token counts in cells.csv.
      3. If a recent run exists (run_000), parse its first snapshot (output00000000_cells.mat) to show counts per cell_type (live + total) and compare.
    """
    import collections, json
    cells_path = CONFIG_FILE.parent / 'cells.csv'
    if not cells_path.exists():
        print('[analyze] cells.csv not found.')
        return
    with open(cells_path,'r') as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if not lines:
        print('[analyze] cells.csv empty.')
        return
    header = lines[0]
    data = lines[1:]
    type_counts = collections.Counter()
    for ln in data:
        parts = ln.split(',')
        if len(parts) < 4:
            continue
        type_counts[parts[3]] += 1
    print(f"[analyze] cells.csv total rows (excluding header): {len(data)}")
    print(f"[analyze] Distinct type labels: {len(type_counts)} -> {dict(list(type_counts.items())[:6])} ...")
    print(f"[analyze] Header columns: {header}")
    if data:
        print(f"[analyze] First 2 lines: {data[:2]}")
        print(f"[analyze] Last 2  lines: {data[-2:]}")

    # Attempt to parse latest run_000 snapshot
    run0 = BASE_OUTPUT_DIR / 'run_000'
    snap_mat = run0 / 'output00000000_cells.mat'
    snap_xml = run0 / 'output00000000.xml'
    if not snap_mat.exists():
        print('[analyze] run_000 snapshot not found; run a test first.')
        return
    if not HAS_SCIPY:
        print('[analyze] scipy missing; cannot open MAT.')
        return
    import numpy as _np, scipy.io as _sio
    try:
        m = _sio.loadmat(snap_mat)
    except Exception as e:  # noqa: BLE001
        print(f"[analyze] MAT load error: {e}")
        return
    arr=None
    for k,v in m.items():
        if k.startswith('__'): continue
        if isinstance(v,_np.ndarray) and v.ndim==2 and v.shape[1]>7:
            arr=v;break
    if arr is None:
        print('[analyze] No cell array found in MAT.')
        return
    ct = arr[:,5].astype(int)
    cyc = arr[:,6] if arr.shape[1]>6 else _np.zeros_like(ct)
    live_mask = (cyc < 100.)
    unique, counts = _np.unique(ct, return_counts=True)
    detail = {}
    for u,c in zip(unique, counts):
        detail[int(u)] = {'total': int(c), 'live': int(((ct==u)&live_mask).sum())}
    print('[analyze] Snapshot0 cell_type distribution (type: {total, live}):')
    print(json.dumps(detail, indent=2))
    # Quick mismatch metric
    print(f"[analyze] Total cells in snapshot0: {int(ct.size)} vs cells.csv lines: {len(data)}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run a parameter sweep for the PhysiCell Recurrence project.")
    parser.add_argument('--test', action='store_true', help="Run in test mode (only one simulation).")
    parser.add_argument('--aggregate', action='store_true', help='After runs, aggregate per-run CSV into population_data.csv (requires pandas).')
    parser.add_argument('--include-dead', action='store_true', help='Include dead cells (cycle_model>=100) in counts.')
    # Removed --threads / --save-interval / --max-time to comply with "only modify specified parameters" requirement.
    parser.add_argument('--use-pymcds', action='store_true', help='Use Studio pyMCDS logic for population counts.')
    parser.add_argument('--microenv', action='store_true', help='When using pyMCDS, also load microenvironment fields (more memory/time).')
    parser.add_argument('--exact-times', action='store_true', help='When using pyMCDS, use each snapshot\'s actual current_time instead of linspace.')
    parser.add_argument('--no-resume', action='store_true', help='Disable resume; run all combinations regardless of existing outputs.')
    parser.add_argument('--redo-failed', action='store_true', help='Re-run runs that have directory but missing population_timeseries.csv.')
    parser.add_argument('--start-run', type=int, default=None, help='Only process runs with run_id >= this value.')
    parser.add_argument('--end-run', type=int, default=None, help='Only process runs with run_id <= this value.')
    parser.add_argument('--analyze-initial', action='store_true', help='Analyze cells.csv vs snapshot0 counts for debugging.')
    parser.add_argument('--cleanup-runs', action='store_true', help='Delete all files in existing run_* folders except population_timeseries.csv and exit.')
    parser.add_argument('--keep-only-population', action='store_true', help='After each successful run, delete all other files in its run_* directory.')
    parser.add_argument('--no-disable-svg', action='store_true', help='Do NOT disable SVG output in temporary settings (default is disabled).')
    args = parser.parse_args()

    if args.cleanup_runs:
        purge_run_directories()
    elif args.analyze_initial:
        analyze_initial_state()
    else:
        main(test_mode=args.test,
             aggregate=args.aggregate,
             include_dead=args.include_dead,
             use_pymcds=args.use_pymcds,
             microenv=args.microenv,
             exact_times=args.exact_times,
             resume=not args.no_resume,
             redo_failed=args.redo_failed,
             start_run=args.start_run,
             end_run=args.end_run,
             keep_only_population=args.keep_only_population,
             disable_svg=not args.no_disable_svg)
