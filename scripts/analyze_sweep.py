#!/usr/bin/env python3
"""Verify all archived sweep scores and regenerate portfolio result figures."""
import csv
import importlib.util
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    project = ROOT / 'models/Recurrence'
    spec = importlib.util.spec_from_file_location('original_scoring', project / 'score_runs.py')
    scoring = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scoring)
    weights = (1, 1, 0.5, 1, 1, 0.5)
    archived = {int(r['run_id']): r for r in csv.DictReader((project / 'scored_runs.csv').open())}
    rows, spans, lengths, max_error = [], set(), set(), 0.0
    for run_id in range(1000):
        path = project / f'sweep_outputs/run_{run_id:03d}/population_timeseries.csv'
        raw = list(csv.DictReader(path.open()))
        times = [float(r['time_min']) for r in raw]
        if any(not math.isfinite(t) for t in times) or any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError(f'Invalid timestamps: {path}')
        for row in raw:
            for key in ('S_aurenus_count', 'engr_S_epidermis_count'):
                value = float(row[key])
                if not math.isfinite(value) or value < 0:
                    raise ValueError(f'Invalid count: {path}')
            if int(row['run_id']) != run_id:
                raise ValueError(f'Run ID mismatch: {path}')
        params = ('degradation_rate', 'aipc_saturation', 'vhh_saturation')
        if any(len({r[k] for r in raw}) != 1 for k in params):
            raise ValueError(f'Parameters vary within a run: {path}')
        row = scoring.score_run(path, weights)
        for key in ('Score', 'P', 'D_A', 'T_A', 'V_A', 'D_E', 'T_E', 'V_E'):
            max_error = max(max_error, abs(float(row[key]) - float(archived[run_id][key])))
        spans.add((times[0], times[-1]))
        lengths.add(len(times))
        rows.append(row)
    if max_error > 1e-12:
        raise ValueError(f'Archived score mismatch: {max_error}')
    rows.sort(key=lambda r: r['Score'], reverse=True)
    best = rows[0]
    report = {
        'runs_verified': len(rows), 'samples_per_run': sorted(lengths),
        'time_spans_minutes': sorted(spans), 'max_archived_score_error': max_error,
        'best_run': best, 'weights': weights,
        'note': 'Verification reuses the original scoring function; it is a consistency check, not independent biological validation.'
    }
    (ROOT / 'results/sweep-verification.json').write_text(json.dumps(report, indent=2) + '\n')
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(projection='3d')
    sc = ax.scatter([float(r['degradation_rate']) for r in rows],
                    [float(r['aipc_saturation']) for r in rows],
                    [float(r['vhh_saturation']) for r in rows],
                    c=[r['Score'] for r in rows], cmap='viridis', s=16, alpha=.75)
    top = rows[:50]
    ax.scatter([float(r['degradation_rate']) for r in top],
               [float(r['aipc_saturation']) for r in top],
               [float(r['vhh_saturation']) for r in top],
               facecolors='none', edgecolors='#dd4a20', s=65, linewidths=1)
    ax.set(xlabel='VHH–cytokine coefficient', ylabel='AIP-C secretion saturation',
           zlabel='VHH secretion saturation', title='1,000 archived parameter combinations\nOrange rings: top 5% by population score')
    fig.colorbar(sc, ax=ax, shrink=.6, pad=.12, label='Score = 1 / (1 + penalty)')
    fig.savefig(ROOT / 'docs/figures/parameter-space.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, row, label in zip(axes, (rows[0], rows[-1]), ('Highest', 'Lowest')):
        rid = int(row['run_id'])
        raw = list(csv.DictReader((project / f'sweep_outputs/run_{rid:03d}/population_timeseries.csv').open()))
        t = [float(r['time_min']) / 60 for r in raw]
        ax.plot(t, [float(r['S_aurenus_count']) for r in raw], color='#6d49a1', label='S. aureus')
        ax.plot(t, [float(r['engr_S_epidermis_count']) for r in raw], color='#19876b', label='Engineered S. epidermidis')
        ax.set(title=f'{label} score: run {rid:03d} ({row["Score"]:.4f})', xlabel='Simulated time (hours)')
        ax.grid(alpha=.2)
    axes[0].set_ylabel('Living agents (microcolony representation)')
    axes[1].legend(loc='best', fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT / 'docs/figures/population-comparison.png', dpi=200)
    plt.close(fig)

    xs = sorted({float(r['degradation_rate']) for r in rows})
    ys = sorted({float(r['vhh_saturation']) for r in rows})
    z = np.array([[np.mean([r['Score'] for r in rows if float(r['degradation_rate']) == x and float(r['vhh_saturation']) == y]) for x in xs] for y in ys])
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(z, origin='lower', aspect='auto', cmap='viridis')
    ax.set_xticks(range(len(xs)), [f'{x:.1f}' for x in xs])
    ax.set_yticks(range(len(ys)), [f'{y:.3f}' for y in ys])
    ax.set(xlabel='VHH–cytokine coefficient', ylabel='VHH secretion saturation',
           title='Mean population score across 10 AIP-C settings')
    fig.colorbar(im, ax=ax, label='Mean score')
    fig.tight_layout()
    fig.savefig(ROOT / 'docs/figures/mean-score.png', dpi=200)
    plt.close(fig)
    print(json.dumps({k: report[k] for k in ('runs_verified', 'samples_per_run', 'time_spans_minutes', 'max_archived_score_error')}))
    print(f'Best run: {best["run_id"]}, score: {best["Score"]:.6f}')


if __name__ == '__main__':
    main()
