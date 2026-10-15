#!/usr/bin/env python3
"""Independently check the archived scores and all 100 mean-score heatmap cells.

Uses the documented equations directly, without importing score_runs.py or
analyze_sweep.py. With --check-image, also checks every heatmap cell's PNG color
against its independently calculated mean; Pillow, NumPy and Matplotlib are needed.
"""
import argparse
import csv
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EPS = 1e-9
COMPONENTS = ('D_A', 'T_A', 'V_A', 'D_E', 'T_E', 'V_E')


def independent_score(raw):
    a = [float(r['S_aurenus_count']) for r in raw]
    e = [float(r['engr_S_epidermis_count']) for r in raw]
    n = len(a)
    if n < 2 or any(not math.isfinite(v) or v < 0 for v in a + e):
        raise ValueError('Invalid archived population series')
    a_scale = a[0] if a[0] > 0 else (math.fsum(a) / n or 1.0)
    e_scale = e[0] if e[0] > 0 else (math.fsum(e) / n or 1.0)
    da = math.fsum(max(v - a[0], 0) / (a_scale + EPS) for v in a) / n
    de = math.fsum(max(e[0] - v, 0) / (e_scale + EPS) for v in e) / n
    ta = max(a[-1] - a[0], 0) / (a_scale + EPS)
    te = max(e[0] - e[-1], 0) / (e_scale + EPS)
    va = math.fsum(max(y - x, 0) for x, y in zip(a, a[1:])) / (a_scale * (n - 1) + EPS)
    ve = math.fsum(max(x - y, 0) for x, y in zip(e, e[1:])) / (e_scale * (n - 1) + EPS)
    penalty = math.fsum((da, ta, .5 * va, de, te, .5 * ve))
    return dict(zip(COMPONENTS, (da, ta, va, de, te, ve)), P=penalty, Score=1 / (1 + penalty))


def verify_png(matrix, image_path):
    import numpy as np
    from PIL import Image
    from matplotlib import colormaps
    from matplotlib.colors import Normalize

    image = np.asarray(Image.open(image_path).convert('RGB'))
    palette = colormaps['viridis'](np.linspace(0, 1, 256), bytes=True)[:, :3]
    def packed(rgb):
        v = rgb.astype(np.int32)
        return (v[..., 0] << 16) | (v[..., 1] << 8) | v[..., 2]
    mask = np.isin(packed(image), packed(palette))
    rows = np.flatnonzero(mask.sum(axis=1) > image.shape[1] * .5)
    if not len(rows):
        raise ValueError('Cannot identify the heatmap color rectangle')
    y0, y1 = int(rows[0]), int(rows[-1]) + 1
    columns = np.flatnonzero(mask[y0:y1].mean(axis=0) > .9)
    groups = np.split(columns, np.flatnonzero(np.diff(columns) > 1) + 1)
    group = max(groups, key=len)
    x0, x1 = int(group[0]), int(group[-1]) + 1
    values = np.asarray(matrix)
    norm = Normalize(vmin=values.min(), vmax=values.max())
    errors = []
    for iy, row in enumerate(matrix):
        for ix, value in enumerate(row):
            # origin='lower': the smallest VHH secretion value is at the bottom.
            x = int(x0 + (ix + .5) * (x1 - x0) / len(row))
            y = int(y0 + (len(matrix) - iy - .5) * (y1 - y0) / len(matrix))
            expected = np.asarray(colormaps['viridis'](norm(value), bytes=True)[:3], dtype=int)
            errors.append(int(np.max(np.abs(image[y, x].astype(int) - expected))))
    max_error = max(errors)
    if max_error > 2:
        raise ValueError(f'Heatmap PNG does not match the independent means: RGB error {max_error}')
    return {'cells_sampled': len(errors), 'max_rgb_channel_error_0_to_255': max_error,
            'heatmap_pixel_bounds_xyxy': [x0, y0, x1, y1], 'smallest_y_at_bottom': True}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check-image', action='store_true')
    ap.add_argument('--output-dir', type=Path, default=ROOT / 'results')
    args = ap.parse_args()
    archive = ROOT / 'models/Recurrence'
    tables = []
    for path in (archive / 'scored_runs.csv', ROOT / 'results/scored_runs.csv'):
        with path.open(newline='') as f:
            rows = list(csv.DictReader(f))
        table = {int(r['run_id']): r for r in rows}
        if len(rows) != 1000 or set(table) != set(range(1000)):
            raise ValueError(f'Incomplete or duplicate score rows: {path}')
        tables.append(table)
    groups = {}
    triplets = set()
    max_error = 0.0
    verified = []
    for rid in range(1000):
        path = archive / f'sweep_outputs/run_{rid:03d}/population_timeseries.csv'
        with path.open(newline='') as f:
            raw = list(csv.DictReader(f))
        if [float(r['time_min']) for r in raw] != list(range(0, 1441, 60)):
            raise ValueError(f'Unexpected timestamps: {path}')
        names = ('degradation_rate', 'aipc_saturation', 'vhh_saturation')
        params = tuple(float(raw[0][name]) for name in names)
        if any(int(r['run_id']) != rid or tuple(float(r[k]) for k in names) != params for r in raw):
            raise ValueError(f'Inconsistent identifiers or parameters: {path}')
        if params in triplets:
            raise ValueError(f'Duplicate parameter combination: {params}')
        triplets.add(params)
        result = independent_score(raw)
        for table in tables:
            if tuple(float(table[rid][k]) for k in names) != params:
                raise ValueError(f'Parameter labels differ from trajectories: run {rid}')
            for key in (*COMPONENTS, 'P', 'Score'):
                max_error = max(max_error, abs(result[key] - float(table[rid][key])))
        groups.setdefault((params[0], params[2]), []).append((params[1], result['Score'], rid))
        verified.append({'run_id': rid, 'coefficient': params[0], 'aipc': params[1],
                         'vhh': params[2], 'Score': result['Score']})
    if max_error > 1e-12:
        raise ValueError(f'Independent score discrepancy: {max_error}')
    xs = sorted({p[0] for p in triplets})
    aips = sorted({p[1] for p in triplets})
    ys = sorted({p[2] for p in triplets})
    if [len(xs), len(aips), len(ys)] != [10, 10, 10] or triplets != set(itertools.product(xs, aips, ys)):
        raise ValueError('The archive is not a complete 10 by 10 by 10 grid')
    matrix, grid_rows = [], []
    for y in ys:
        matrix_row = []
        for x in xs:
            items = sorted(groups[x, y])
            if [v[0] for v in items] != aips:
                raise ValueError('Missing or duplicated AIP-C setting in a mean-score cell')
            mean = math.fsum(v[1] for v in items) / len(items)
            matrix_row.append(mean)
            grid_rows.append({'degradation_rate': x, 'vhh_saturation': y, 'n_aipc_settings': len(items),
                              'mean_score': mean, 'min_score': min(v[1] for v in items),
                              'max_score': max(v[1] for v in items),
                              'run_ids': ';'.join(str(v[2]) for v in items)})
        matrix.append(matrix_row)
    report = {'verification_added': '2026-10-15', 'independent_of_original_scorer': True,
              'runs_recomputed': len(verified), 'max_score_or_component_difference': max_error,
              'grid_shape': [len(ys), len(xs)], 'aipc_settings_per_cell': 10,
              'x_values': xs, 'y_values': ys, 'aipc_values': aips,
              'matrix_rows_in_ascending_y': matrix,
              'highest_mean_cell': max(grid_rows, key=lambda r: r['mean_score']),
              'lowest_mean_cell': min(grid_rows, key=lambda r: r['mean_score']),
              'note': 'Arithmetic mean of 10 population scores per cell; each parameter combination has one archived trajectory. This does not validate biological outcomes or stochastic robustness.'}
    if args.check_image:
        report['existing_png_check'] = verify_png(matrix, ROOT / 'docs/figures/mean-score.png')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / 'mean-score-grid.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(grid_rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(grid_rows)
    (args.output_dir / 'mean-score-verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['runs_recomputed', 'max_score_or_component_difference',
                                           'highest_mean_cell', 'lowest_mean_cell']}, indent=2))
    if args.check_image:
        print(json.dumps(report['existing_png_check']))


if __name__ == '__main__':
    main()
