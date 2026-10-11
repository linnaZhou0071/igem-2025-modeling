#!/usr/bin/env python3
"""
Score all sweep runs (0-999 by default) using the multi-component penalty -> score formula:

Components (directional, only penalize undesired direction):
  D_A, D_E : average baseline deviation (area-style)
  T_A, T_E : net trend (start->end)
  V_A, V_E : cumulative directional step violations

Penalty: P = w_A^D D_A + w_A^T T_A + w_A^V V_A + w_E^D D_E + w_E^T T_E + w_E^V V_E
Score  : 1 / (1 + P)   (gamma=0, no ratio reward R in this initial implementation)

Outputs a CSV (default: scored_runs.csv) with detailed components for inspection.

Usage examples:
  python score_runs.py                       # score runs 0..999 using defaults
  python score_runs.py --min-run 0 --max-run 199 --output scored_0_199.csv
  python score_runs.py --verify-run 4        # print detailed breakdown for run_004
  python score_runs.py --weights 1,1,0.5,1,1,0.5  # (w_A^D,w_A^T,w_A^V,w_E^D,w_E^T,w_E^V)

Verification approach suggestions (also printed when using --verify-run):
  1. Manually recompute one row with a quick inline Python snippet (provided).
  2. Check that V_A * (A0*(N-1)) approximates sum of all positive forward increments.
  3. Confirm D_A equals mean( max(0,(A_i-A0)/A0) ). Same for E with direction inverted.
  4. Ensure T_A matches max(0,(A_end-A0)/A0) and T_E matches max(0,(E0-E_end)/E0).

Edge handling:
  - If initial count is 0, we fall back to using mean(series) for denominators to avoid zero-division inflation.
  - Missing run folders or CSV files are skipped with a warning.

Author: automated assistant
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple

EPS = 1e-9

EXPECTED_FILENAME = "population_timeseries.csv"
# Column names expected in the CSV
COL_TIME = "time_min"
COL_A = "S_aurenus_count"
COL_E = "engr_S_epidermis_count"
COL_RUN_ID = "run_id"  # we still trust file path for run id
COL_DEGR = "degradation_rate"
COL_AIPC = "aipc_saturation"
COL_VHH = "vhh_saturation"

class RunScoreError(Exception):
    pass

def safe_int_or_float(v: str) -> float:
    try:
        if v.strip() == '':
            return 0.0
        if '.' in v or 'e' in v.lower():
            return float(v)
        return float(int(v))
    except ValueError:
        return float(v)

def compute_components(A: List[float], E: List[float]) -> Dict[str, float]:
    if not A or not E:
        raise RunScoreError("Empty series for A or E")
    if len(A) != len(E):
        raise RunScoreError("Mismatched lengths for A and E series")

    N = len(A)
    A0, E0 = A[0], E[0]
    A_end, E_end = A[-1], E[-1]

    # Fallback denominators if initial is zero.
    base_A = A0 if A0 > 0 else (sum(A)/N if sum(A) > 0 else 1.0)
    base_E = E0 if E0 > 0 else (sum(E)/N if sum(E) > 0 else 1.0)

    # 1. Baseline deviation (directional)
    D_A = sum(max(0.0, (a - A0) / (A0 + EPS)) for a in A) / N if A0 > 0 else sum(max(0.0, (a - A0) / (base_A + EPS)) for a in A) / N
    D_E = sum(max(0.0, (E0 - e) / (E0 + EPS)) for e in E) / N if E0 > 0 else sum(max(0.0, (E0 - e) / (base_E + EPS)) for e in E) / N

    # 2. Net trend (start->end, directional)
    T_A = max(0.0, (A_end - A0) / (A0 + EPS)) if A0 > 0 else max(0.0, (A_end - A0) / (base_A + EPS))
    T_E = max(0.0, (E0 - E_end) / (E0 + EPS)) if E0 > 0 else max(0.0, (E0 - E_end) / (base_E + EPS))

    # 3. Directional single-step violations
    pos_steps_A = 0.0
    neg_steps_E = 0.0
    for i in range(1, N):
        dA = A[i] - A[i-1]
        if dA > 0:
            pos_steps_A += dA
        dE = E[i-1] - E[i]
        if dE > 0:
            neg_steps_E += dE
    denom_A_steps = (A0 if A0 > 0 else base_A) * (N - 1) + EPS
    denom_E_steps = (E0 if E0 > 0 else base_E) * (N - 1) + EPS
    V_A = pos_steps_A / denom_A_steps
    V_E = neg_steps_E / denom_E_steps

    return dict(
        N=N,
        A0=A0, E0=E0,
        A_end=A_end, E_end=E_end,
        D_A=D_A, D_E=D_E,
        T_A=T_A, T_E=T_E,
        V_A=V_A, V_E=V_E,
        pos_steps_A=pos_steps_A,
        neg_steps_E=neg_steps_E,
    )

def score_run(csv_path: Path, weights: Tuple[float,float,float,float,float,float]):
    w_A_D, w_A_T, w_A_V, w_E_D, w_E_T, w_E_V = weights
    if not csv_path.is_file():
        raise RunScoreError(f"Missing file: {csv_path}")

    with csv_path.open('r', newline='') as f:
        reader = csv.DictReader(f)
        A: List[float] = []
        E: List[float] = []
        params_last = {}  # keep last row's parameters (they should be identical across rows)
        for row in reader:
            try:
                A.append(safe_int_or_float(row[COL_A]))
                E.append(safe_int_or_float(row[COL_E]))
                params_last = {
                    'run_id': row.get(COL_RUN_ID, ''),
                    'degradation_rate': row.get(COL_DEGR, ''),
                    'aipc_saturation': row.get(COL_AIPC, ''),
                    'vhh_saturation': row.get(COL_VHH, ''),
                }
            except KeyError as ke:
                raise RunScoreError(f"Missing expected column {ke} in {csv_path}") from ke

    comps = compute_components(A, E)
    P = (
        w_A_D * comps['D_A'] +
        w_A_T * comps['T_A'] +
        w_A_V * comps['V_A'] +
        w_E_D * comps['D_E'] +
        w_E_T * comps['T_E'] +
        w_E_V * comps['V_E']
    )
    score = 1.0 / (1.0 + P)  # gamma=0, no reward term

    return {
        **params_last,
        **comps,
        'P': P,
        'Score': score,
    }

def parse_weights(s: str) -> Tuple[float,float,float,float,float,float]:
    parts = [p.strip() for p in s.split(',') if p.strip()]
    if len(parts) != 6:
        raise argparse.ArgumentTypeError("--weights requires 6 comma-separated numbers")
    try:
        vals = tuple(float(p) for p in parts)  # type: ignore
    except ValueError as e:
        raise argparse.ArgumentTypeError("Invalid float in --weights") from e
    return vals  # type: ignore

def main():
    ap = argparse.ArgumentParser(description="Score parameter sweep runs using multi-component penalty -> score formula.")
    ap.add_argument('--sweep-root', type=Path, default=Path(__file__).parent / 'sweep_outputs', help='Root folder containing run_XXX directories')
    ap.add_argument('--min-run', type=int, default=0)
    ap.add_argument('--max-run', type=int, default=999)
    ap.add_argument('--output', type=Path, default=Path(__file__).parent / 'scored_runs.csv')
    ap.add_argument('--weights', type=parse_weights, default=(1.0,1.0,0.5, 1.0,1.0,0.5), help='Comma list: w_A^D,w_A^T,w_A^V,w_E^D,w_E^T,w_E^V')
    ap.add_argument('--verify-run', type=int, default=None, help='Print detailed breakdown for a single run id (after scoring all)')
    ap.add_argument('--skip-missing', action='store_true', help='Silently skip runs missing CSV file')
    ap.add_argument('--verbose', '-v', action='store_true')
    args = ap.parse_args()

    sweep_root: Path = args.sweep_root
    if not sweep_root.exists():
        print(f"[ERROR] Sweep root not found: {sweep_root}", file=sys.stderr)
        sys.exit(1)

    rows_out = []
    failures = 0
    total = 0
    for rid in range(args.min_run, args.max_run + 1):
        run_dir = sweep_root / f"run_{rid:03d}"
        csv_file = run_dir / EXPECTED_FILENAME
        total += 1
        try:
            result = score_run(csv_file, args.weights)
            # ensure run_id if missing
            if not result.get('run_id'):
                result['run_id'] = str(rid)
            result['run_dir'] = run_dir.name
            rows_out.append(result)
            if args.verbose and rid % 50 == 0:
                print(f"[info] processed run {rid:03d}")
        except Exception as e:
            failures += 1
            msg = f"[warn] run {rid:03d} skipped: {e}"
            if args.skip-missing:
                if args.verbose:
                    print(msg)
            else:
                print(msg, file=sys.stderr)

    # Sort by Score descending
    rows_out.sort(key=lambda d: d['Score'], reverse=True)

    # Write output CSV
    if rows_out:
        fieldnames = [
            'rank','Score','P',
            'run_id','run_dir',
            'degradation_rate','aipc_saturation','vhh_saturation',
            'N','A0','A_end','E0','E_end',
            'D_A','T_A','V_A','D_E','T_E','V_E',
            'pos_steps_A','neg_steps_E'
        ]
        with args.output.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for rank, row in enumerate(rows_out, start=1):
                out_row = {k: row.get(k, '') for k in fieldnames}
                out_row['rank'] = rank
                writer.writerow(out_row)
        print(f"[done] Wrote {len(rows_out)} scored runs to {args.output}")
        if failures:
            print(f"[info] Failures/skipped: {failures}/{total}")
    else:
        print("[warn] No runs scored; check paths or patterns")

    # Verification breakdown if requested
    if args.verify_run is not None:
        rid = args.verify_run
        match = None
        for r in rows_out:
            # run_id might be string, ensure equality
            if str(r.get('run_id')) == str(rid):
                match = r
                break
        if not match:
            print(f"[verify] Run {rid} not found in scored rows (maybe missing file or outside range)")
        else:
            print("\n[verify] Detailed components for run", rid)
            for k in ['Score','P','D_A','T_A','V_A','D_E','T_E','V_E','A0','A_end','E0','E_end','pos_steps_A','neg_steps_E','N']:
                print(f"  {k}: {match[k]}")
            print("[verify] Manual one-liner check template (adjust run id/path if needed):\n")
            print(f"python - <<'PY'\nimport csv; import math; f='{(sweep_root / f'run_{rid:03d}' / EXPECTED_FILENAME).as_posix()}';\nrows=list(csv.DictReader(open(f))); A=[float(r['{COL_A}']) for r in rows]; E=[float(r['{COL_E}']) for r in rows];\nA0=A[0]; E0=E[0]; eps=1e-9; N=len(A)\nD_A=sum(max(0,(a-A0)/(A0+eps)) for a in A)/N\nD_E=sum(max(0,(E0-e)/(E0+eps)) for e in E)/N\nT_A=max(0,(A[-1]-A0)/(A0+eps)); T_E=max(0,(E0-E[-1])/(E0+eps))\npos_steps=sum(max(0,A[i]-A[i-1]) for i in range(1,N))\nneg_steps=sum(max(0,E[i-1]-E[i]) for i in range(1,N))\nV_A=pos_steps/(A0*(N-1)+eps); V_E=neg_steps/(E0*(N-1)+eps)\nP = (1*D_A + 1*T_A + 0.5*V_A + 1*D_E + 1*T_E + 0.5*V_E)\nScore=1/(1+P)\nprint('D_A',D_A,'T_A',T_A,'V_A',V_A,'D_E',D_E,'T_E',T_E,'V_E',V_E,'P',P,'Score',Score)\nPY")

if __name__ == '__main__':
    main()
