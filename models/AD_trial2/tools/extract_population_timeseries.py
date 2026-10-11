#!/usr/bin/env python3
"""extract_population_timeseries.py

Extract per-cell-type population counts vs time from an existing PhysiCell
project output folder (no re-simulation). Designed for the AD_trial2 project.

Logic overview:
1. Discover output snapshot basenames: outputXXXXXXXX.xml (sorted ascending).
2. For each snapshot, parse time from <current_time> in the XML (fallback: index * --fallback-interval if provided).
3. Load the corresponding outputXXXXXXXX_cells.mat file for cell data.
4. Detect the cell type row heuristically and count cells for each defined cell type ID.
5. (Optional) Exclude dead cells by detecting a cycle / death state row (values >=100 treated as dead), unless --include-dead.
6. Write a CSV in a format similar to population_timeseries.csv used in parameter sweeps.

Header format:
time_min,<typeName1>_count,<typeName2>_count,...,run_id[,<extra param columns>]

Cell type names are parsed from the first snapshot XML (<cell_definitions><cell_definition name="..." ID="...">).
If parsing fails, generic names type<ID> are used.

Heuristics for detecting the cell_type row in the MAT file:
 - Load candidate 2D arrays where min(dim) < max(dim); orient as (features, n_cells).
 - Among rows, find integer-like rows. Pick the one whose unique integer set overlaps most with known cell definition IDs.
 - If no mapping available, pick integer-like row with highest unique count (<= 32 distinct values).

Cycle/live detection (for excluding dead cells unless --include-dead):
 - Search for another integer-like row (not the type row) containing some values >=100; treat values >=100 as dead.
 - If not found, assume all cells are live.

Usage examples:
  python tools/extract_population_timeseries.py --output-dir output \
         --out-csv population_timeseries.csv

  python tools/extract_population_timeseries.py --output-dir output \
         --out-csv population_timeseries.csv --param degradation_rate=0.1 \
         --param aipc_saturation=0.5 --param vhh_saturation=0.01

  python tools/extract_population_timeseries.py --output-dir output \
         --out-csv pop_dead_included.csv --include-dead

"""
from __future__ import annotations

import argparse
import os
import re
import sys
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
try:
    import scipy.io as sio
except ImportError as e:  # noqa: BLE001
    print('[error] scipy is required for this script (missing).')
    raise
import xml.etree.ElementTree as ET


def parse_args():
    p = argparse.ArgumentParser(description='Extract per-cell-type population time series from existing PhysiCell outputs.')
    p.add_argument('--output-dir', default='output', help='Directory containing output*.xml and *_cells.mat files (default: output)')
    p.add_argument('--out-csv', default='population_timeseries.csv', help='Result CSV filename (default: population_timeseries.csv)')
    p.add_argument('--include-dead', action='store_true', help='Include dead cells (values with cycle_model >= 100) in counts')
    p.add_argument('--run-id', type=int, default=0, help='Run identifier to write in CSV (default: 0)')
    p.add_argument('--param', action='append', default=[], metavar='KEY=VAL', help='Extra constant parameter columns to append (e.g. --param degradation_rate=0.1)')
    p.add_argument('--fallback-interval', type=float, default=None, help='If a time cannot be parsed from XML, use index * this value (minutes).')
    p.add_argument('--id-filter', nargs='*', help='Optional subset of cell type names (exact) to include; default: all.')
    p.add_argument('--id-alias', action='append', default=[], metavar='ID=Name', help='Override or supply a name for a numeric cell type ID (e.g. --id-alias 1=macrophage). Multiple allowed.')
    p.add_argument('--list-only', action='store_true', help='Only list discovered / aliased cell type mapping then exit (no CSV).')
    p.add_argument('--settings', default=None, help='Optional path to settings.xml to supplement/override cell definition names (useful when snapshot XML lacks full <cell_definitions>).')
    p.add_argument('--verbose', action='store_true', help='Verbose logging.')
    return p.parse_args()


def discover_snapshots(output_dir: Path) -> List[str]:
    snaps = sorted([f for f in os.listdir(output_dir) if f.startswith('output') and f.endswith('.xml')])
    return snaps


def parse_time_from_xml(xml_path: Path) -> Optional[float]:
    pattern = re.compile(r'<current_time[^>]*>([-+0-9.eE]+)')
    try:
        with open(xml_path, 'r') as f:
            for line in f:
                m = pattern.search(line)
                if m:
                    try:
                        return float(m.group(1))
                    except Exception:
                        return None
    except OSError:
        return None
    return None


def parse_cell_definitions(xml_path: Path) -> Dict[int, str]:
    """Parse cell definitions robustly.

    Handles variations: attributes may be 'ID', 'id', 'type', 'typeID'. Name may be in 'name' or 'type'.
    Returns mapping of numeric ID -> name.
    """
    mapping: Dict[int, str] = {}
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        for cdef in root.findall('.//cell_definitions/cell_definition'):
            # Candidate id attributes
            cid_txt = (cdef.get('ID') or cdef.get('id') or cdef.get('type') or cdef.get('typeID'))
            if cid_txt is None:
                continue
            try:
                cid = int(cid_txt)
            except ValueError:
                continue
            name = cdef.get('name') or cdef.get('type') or f'type{cid}'
            mapping[cid] = name
    except Exception:
        return mapping
    return mapping


def locate_settings_file(start_dir: Path) -> Optional[Path]:
    """Search upward from start_dir for a config/settings.xml file if user didn't supply one."""
    candidates = []
    for up in [start_dir] + list(start_dir.parents):
        cfg = up / 'config' / 'settings.xml'
        if cfg.exists():
            candidates.append(cfg)
    # prefer the closest (smallest depth difference)
    if candidates:
        return candidates[0]
    return None


def sanitize_name(name: str) -> str:
    import re
    s = re.sub(r'[^0-9A-Za-z]+', '_', name.strip())
    s = re.sub(r'_+', '_', s).strip('_')
    return s if s else 'unnamed'


def load_cells_matrix(mat_path: Path) -> Optional[np.ndarray]:
    try:
        m = sio.loadmat(mat_path)
    except Exception as e:  # noqa: BLE001
        print(f'[warn] could not load {mat_path.name}: {e}')
        return None
    best = None
    for k, v in m.items():
        if k.startswith('__'):
            continue
        if isinstance(v, np.ndarray) and v.ndim == 2:
            r, c = v.shape
            if min(r, c) < max(r, c):  # rectangular enough
                # orient as (features, n_cells)
                if r <= c:
                    cand = v
                else:
                    cand = v.T
                # heuristic: need at least a few columns
                if cand.shape[1] >= 10:
                    best = cand
                    break
    return best


def detect_type_and_cycle_rows(cells: np.ndarray) -> Tuple[int, Optional[int]]:
    """Replicate the detection heuristic used in parameter_sweep.extract_population_data.

    We scan rows (features) of the (features, n_cells) matrix. For each row:
      - Compute unique values u.
      - If 2 <= len(u) <= 10 and all integer-like:
          * If set(u) includes {0,5}, mark as potential cell_type row (ct_row).
          * Else if cycle_row not set yet, mark as potential cycle row (cyc_row).
    Fallbacks: ct_row -> 5 if not found, cyc_row -> 4 if not found (mirrors original script).
    """
    f, n = cells.shape
    ct_row = None
    cyc_row = None
    for ri in range(f):
        row = cells[ri, :]
        if not np.allclose(row, np.round(row), atol=1e-6):
            continue
        u = np.unique(row.astype(int))
        if 2 <= len(u) <= 10:
            if set(u).issuperset({0, 5}):
                ct_row = ri
            elif cyc_row is None:
                cyc_row = ri
    if ct_row is None:
        ct_row = 5 if f > 5 else 0
    if cyc_row is None:
        cyc_row = 4 if f > 4 else None
    return ct_row, cyc_row


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    if not out_dir.exists():
        print(f'[error] output directory {out_dir} does not exist.')
        sys.exit(1)

    snapshots = discover_snapshots(out_dir)
    if not snapshots:
        print(f'[error] no output*.xml files found in {out_dir}')
        sys.exit(1)

    # Parse cell definitions from first snapshot
    first_xml = out_dir / snapshots[0]
    id_to_name = parse_cell_definitions(first_xml)
    if args.verbose:
        print(f'[info] Parsed {len(id_to_name)} cell definitions from first snapshot XML.')

    # If user supplied settings.xml, merge/override names
    if args.settings:
        settings_path = Path(args.settings)
        if settings_path.exists():
            extra_defs = parse_cell_definitions(settings_path)
            if args.verbose:
                print(f'[info] Parsed {len(extra_defs)} cell definitions from provided settings.xml.')
            for cid, nm in extra_defs.items():
                # do not overwrite explicit alias later; but replace placeholder typeX
                if (cid not in id_to_name) or id_to_name[cid].startswith('type'):
                    id_to_name[cid] = nm
        else:
            print(f'[warn] --settings path not found: {settings_path}')
    else:
        # Auto-locate settings if we got suspiciously few definitions (<2 or missing sequential IDs)
        if len(id_to_name) < 3 or not all(i in id_to_name for i in range(max(id_to_name.keys())+1)):
            located = locate_settings_file(out_dir)
            if located and located != first_xml:
                extra_defs = parse_cell_definitions(located)
                if args.verbose:
                    print(f'[info] Supplemented {len(extra_defs)} definitions from {located}')
                for cid, nm in extra_defs.items():
                    if (cid not in id_to_name) or id_to_name[cid].startswith('type'):
                        id_to_name[cid] = nm
    known_ids = sorted(id_to_name.keys())

    # Apply manual aliases (override or add)
    for alias in args.id_alias:
        if '=' in alias:
            k,v = alias.split('=',1)
            try:
                cid = int(k)
                id_to_name[cid] = v.strip()
                if cid not in known_ids:
                    known_ids.append(cid)
            except ValueError:
                print(f'[warn] ignoring id-alias (non-integer id): {alias}')
        else:
            print(f'[warn] ignoring malformed --id-alias {alias} (expected ID=Name)')
    known_ids = sorted(set(known_ids))

    if args.list_only:
        print('Discovered / aliased cell type mapping:')
        for cid in known_ids:
            print(f'  {cid}: {id_to_name.get(cid, f"type{cid}")}')
        return
    if not known_ids and args.verbose:
        print('[info] no cell definitions parsed; will derive names from encountered IDs.')

    # Parameter columns
    param_kv: Dict[str, str] = {}
    for tok in args.param:
        if '=' in tok:
            k, v = tok.split('=', 1)
            param_kv[k.strip()] = v.strip()
        else:
            print(f'[warn] ignoring malformed --param {tok} (expected key=value)')

    # Data accumulators
    time_points: List[float] = []
    counts_per_type: Dict[int, List[int]] = {}  # id -> list aligned with time_points
    # Pre-seed counts_per_type with all defined IDs so even absent types get zero columns
    for _cid in known_ids:
        counts_per_type[_cid] = []

    # Prepare cycle/live detection once we load first MAT
    cell_type_row: Optional[int] = None
    cycle_row: Optional[int] = None

    for idx, xml_name in enumerate(snapshots):
        base = xml_name.replace('.xml', '')
        mat_path = out_dir / f'{base}_cells.mat'
        xml_path = out_dir / xml_name
        if not mat_path.exists():
            if args.verbose:
                print(f'[skip] {mat_path.name} missing.')
            continue
        # time
        t = parse_time_from_xml(xml_path)
        if t is None and args.fallback_interval is not None:
            # derive index from base (outputXXXXXXXX)
            try:
                serial = int(base.replace('output', ''))
                t = serial * args.fallback_interval
            except Exception:
                t = 0.0
        elif t is None:
            t = 0.0

        cells = load_cells_matrix(mat_path)
        if cells is None:
            continue
        if cell_type_row is None:
            cell_type_row, cycle_row = detect_type_and_cycle_rows(cells)
            if args.verbose:
                print(f'[detect] cell_type_row={cell_type_row}; cycle_row={cycle_row}')

        ct_vals = cells[cell_type_row, :].astype(int)
        if cycle_row is not None and not args.include_dead:
            cyc = cells[cycle_row, :]
            live_mask = (cyc < 100.)
        else:
            live_mask = np.ones_like(ct_vals, dtype=bool)

        # Expand ID->name mapping with any new IDs appearing dynamically (very rare)
        new_ids = np.unique(ct_vals)
        added_dynamic = False
        for cid in new_ids:
            if cid not in id_to_name:
                id_to_name[cid] = f'type{cid}'
                counts_per_type[cid] = [0]*len(time_points)  # back-fill zeros for previous snapshots
                added_dynamic = True
        # Canonical names for 0 / 5 if not aliased
        if 0 in id_to_name and id_to_name.get(0,'').startswith('type'):
            id_to_name[0] = 'S_aurenus'
        if 5 in id_to_name and id_to_name.get(5,'').startswith('type'):
            id_to_name[5] = 'engr_S_epidermis'

        # Snapshot index
        snapshot_index = len(time_points)

        # If a new type appears mid-run, pad its previous snapshots with zeros
        for cid in id_to_name.keys():
            if cid not in counts_per_type:
                counts_per_type[cid] = [0] * snapshot_index

        # Count live cells for each type in this snapshot
        for cid in counts_per_type.keys():
            mask = (ct_vals == cid) & live_mask
            counts_per_type[cid].append(int(mask.sum()))

        time_points.append(t)

    if not time_points:
        print('[error] no valid snapshots processed.')
        sys.exit(1)

    # Normalize all count lists to equal length
    T = len(time_points)
    for cid, lst in counts_per_type.items():
        if len(lst) < T:
            lst.extend([0] * (T - len(lst)))

    # Decide header order (respect given id_filter order if provided)
    if args.id_filter:
        # Keep only those present; preserve filter order
        ordered_ids = []
        name_to_id = {nm: cid for cid, nm in id_to_name.items()}
        for nm in args.id_filter:
            if nm in name_to_id:
                ordered_ids.append(name_to_id[nm])
            else:
                print(f'[warn] id-filter name not found: {nm}')
    else:
        ordered_ids = sorted(id_to_name.keys())

    # Build header
    header_names = [f'{sanitize_name(id_to_name[cid])}_count' for cid in ordered_ids]
    param_cols = list(param_kv.keys())
    out_path = Path(args.out_csv)
    with open(out_path, 'w') as f:
        f.write('time_min,')
        f.write(','.join(header_names))
        f.write(',run_id')
        if param_cols:
            for c in param_cols:
                f.write(f',{c}')
        f.write('\n')
        for ti, t in enumerate(time_points):
            row_counts = [str(counts_per_type[cid][ti]) for cid in ordered_ids]
            line = f"{t},{','.join(row_counts)},{args.run_id}"
            if param_cols:
                for c in param_cols:
                    line += f',{param_kv[c]}'
            f.write(line + '\n')
    print(f'[done] Wrote {out_path} with {len(time_points)} time points and {len(ordered_ids)} cell type columns.')

    if args.verbose:
        print('Cell type mapping:')
        for cid in ordered_ids:
            print(f'  ID {cid} -> {id_to_name[cid]}')


if __name__ == '__main__':
    main()
