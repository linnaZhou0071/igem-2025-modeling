#!/usr/bin/env python3
"""
Generate modified cells.csv with a single centered sebaceous gland based on prior Variant D logic.
Changes requested:
 1. Only one gland, centered in the xy domain (computed from keratinocyte bounding box midpoint).
 2. Implement a real spherical cap (multiple z levels) instead of a single flat bottom disc.
 3. Reduce bottom radius to 100 µm (top radius remains 50 µm; neck length 40 µm; power ramp to 100 µm at depth 150 µm).

Assumptions (can be adjusted easily):
 - TotalDepth kept at 150 µm (unchanged) for the cylindrical + taper section.
 - Spherical cap extends an additional CapHeight (default 20 µm) deeper than the nominal bottom to provide curvature closure.
   (If you prefer not to extend deeper, set CapHeight=0 or adjust.)
 - Cap radius profile: hemisphere fragment of radius BottomRadius: r(z_extra) = sqrt(BottomRadius^2 - z_extra^2) for z_extra in [0, CapHeight].
   We truncate at CapHeight << BottomRadius so closure reaches near 0 radius without excessive depth increase.
 - Removal logic purges keratinocytes inside the gland volume for depths within TotalDepth (not including cap extension) to avoid over-removal deeper.

Outputs:
 - cells_with_single_gland_D.csv adjacent to original cells.csv
"""
from __future__ import annotations
import csv, math, os, sys
from statistics import mean

BASE_CELLS = os.path.join(os.path.dirname(__file__), "..", "config", "cells.csv")
OUT_SINGLE = os.path.join(os.path.dirname(BASE_CELLS), "cells_with_single_gland_D.csv")

TopRadius = 50.0
BottomRadius = 80.0   # new reduced bottom radius (was 100)
TotalDepth = 150.0
NeckDepth = 40.0
PowerExp = 1.2
LayerStep = 10.0
ShellMargin = 4.0
CircumSpacing = 15.0
BendAmp = 0.0  # No bending for single centered gland unless desired
# CapHeight request: reduce extension depth; with BottomRadius=80 the maximum hemispherical
# cap height is 80 (cap height cannot exceed sphere radius). Using full hemisphere again; if
# you want a truncated cap (keeping a non-point bottom) set CapHeight < BottomRadius.
CapHeight = 80.0
CapStep = 5.0             # finer vertical spacing for smoother cap & guaranteed closure
BottomRefineThreshold = 20.0  # start refining steps when remaining depth < this (µm)
BottomRefineStep = 1.0        # finer step near bottom (option B)
SmallRadiusClose = 2.0        # tighter closure threshold (option A)

# ----------------- IO helpers -----------------

def read_cells(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]

def parse_row(row):
    x = float(row[0]); y = float(row[1]); z = float(row[2]); ctype = row[3].strip()
    return x,y,z,ctype

# ----------------- Geometry detection -----------------

def find_third_layer_z(data):
    z_vals = []
    for r in data:
        try:
            _,_,z,ctype = parse_row(r)
        except Exception:
            continue
        if ctype == 'keratinocyte':
            z_vals.append(round(z))
    uniq = sorted(set(z_vals), reverse=True)
    if len(uniq) < 3:
        raise RuntimeError(f"Not enough keratinocyte z-layers: {uniq}")
    return float(uniq[2])


def compute_center(data):
    xs, ys = [], []
    for r in data:
        try:
            x,y,_,ctype = parse_row(r)
        except Exception:
            continue
        if ctype == 'keratinocyte':
            xs.append(x); ys.append(y)
    if not xs:
        raise RuntimeError("No keratinocytes found to compute center.")
    # center of bounding box, not average (less sensitive to density variations)
    cx = 0.5*(min(xs) + max(xs))
    cy = 0.5*(min(ys) + max(ys))
    return cx, cy

# ----------------- Radius profile -----------------

def radius_profile(depth):
    if depth <= NeckDepth:
        return TopRadius
    if depth > TotalDepth:
        return 0.0
    frac = (depth - NeckDepth) / (TotalDepth - NeckDepth)
    return TopRadius + (BottomRadius - TopRadius) * (frac ** PowerExp)

# ----------------- Volume tests -----------------

def extended_radius(depth):
    """Radius including cap region for interior clearing."""
    if depth < 0:
        return 0.0
    if depth <= TotalDepth:
        return radius_profile(depth)
    cap_limit = TotalDepth + CapHeight
    if depth <= cap_limit:
        extra = depth - TotalDepth
        # hemispherical cap; protect against floating error
        return math.sqrt(max(0.0, BottomRadius**2 - extra**2))
    return 0.0

def inside_gland_volume(x,y,z,z_open,center):
    depth = z_open - z
    # allow purge through cap region
    if depth < 0 or depth > (TotalDepth + CapHeight):
        return False
    cx,cy = center
    r = extended_radius(depth)
    r_allow = r - ShellMargin
    if r_allow <= 0:
        return False
    if (x-cx)**2 + (y-cy)**2 < r_allow**2:
        return True
    return False

# ----------------- Shell construction -----------------

def build_shell(z_open, center):
    cx,cy = center
    shell_pts = []
    depths = [d for d in [i*LayerStep for i in range(int(TotalDepth/LayerStep)+1)] if d <= TotalDepth]
    for d in depths:
        r = radius_profile(d)
        if r <= 0: continue
        n_theta = max(6, round(2*math.pi*r / CircumSpacing))
        z = z_open - d
        for k in range(n_theta):
            theta = 2*math.pi*k / n_theta
            x = cx + r*math.cos(theta)
            y = cy + r*math.sin(theta)
            shell_pts.append((x,y,z))
    # Spherical cap closure (extends deeper by CapHeight)
    if CapHeight > 0:
        # dynamic stepping toward bottom for denser closure
        z_extra = 0.0
        while z_extra < CapHeight - 1e-6:
            remaining = CapHeight - z_extra
            step = CapStep if remaining > BottomRefineThreshold else BottomRefineStep
            # ensure we land exactly on CapHeight
            if z_extra + step > CapHeight:
                step = CapHeight - z_extra
            z_extra += step
            r_cap = math.sqrt(max(0.0, BottomRadius**2 - z_extra**2))
            z = (z_open - TotalDepth) - z_extra
            if r_cap <= SmallRadiusClose:
                # create a small ring (if r_cap > ~0) at same depth, plus centerpoint, then stop
                if r_cap > 0.5:
                    n_theta = max(6, round(2*math.pi*max(r_cap,1.0) / CircumSpacing))
                    for k in range(n_theta):
                        theta = 2*math.pi*k / n_theta
                        x = cx + r_cap*math.cos(theta)
                        y = cy + r_cap*math.sin(theta)
                        shell_pts.append((x,y,z))
                # center point at same z for sealed closure
                shell_pts.append((cx,cy,z))
                break
            n_theta = max(6, round(2*math.pi*r_cap / CircumSpacing))
            for k in range(n_theta):
                theta = 2*math.pi*k / n_theta
                x = cx + r_cap*math.cos(theta)
                y = cy + r_cap*math.sin(theta)
                shell_pts.append((x,y,z))
    return shell_pts

# ----------------- Processing -----------------

def process(header, data, z_open, center):
    new_rows = []
    cut_removed = 0
    vol_removed = 0
    cx,cy = center
    top_r2 = TopRadius**2
    for row in data:
        try:
            x,y,z,ctype = parse_row(row)
        except Exception:
            continue
        if ctype == 'keratinocyte':
            # initial top cylinder cut (only one center)
            if (x-cx)**2 + (y-cy)**2 < top_r2:
                cut_removed += 1
                continue
            if inside_gland_volume(x,y,z,z_open,center):
                vol_removed += 1
                continue
        new_rows.append(row)
    shell = build_shell(z_open, center)
    for (x,y,z) in shell:
        new_rows.append([f"{x:.3f}", f"{y:.3f}", f"{z:.3f}", 'keratinocyte', '', '', '', ''])
    # --- Inject sebocytes at bottom interior ---
    sebocytes_added = 0
    # Estimate geometric bottom (even if closure earlier):
    z_bottom = z_open - (TotalDepth + CapHeight)
    # Place them slightly above true bottom to avoid merging with closure center point
    z_seb = z_bottom + 4.0  # 4 µm above bottom
    r_seb = 6.0             # radial distance from center
    seb_positions = [
        (cx + r_seb, cy, z_seb),
        (cx - r_seb, cy, z_seb),
        (cx, cy + r_seb, z_seb),
        (cx, cy - r_seb, z_seb)
    ]
    for (x,y,z) in seb_positions:
        new_rows.append([f"{x:.3f}", f"{y:.3f}", f"{z:.3f}", 'sebocyte', '', '', '', ''])
        sebocytes_added += 1
    stats = {
        'cut_removed': cut_removed,
        'vol_removed': vol_removed,
        'shell_added': len(shell),
        'sebocytes_added': sebocytes_added,
        'final_total': len(new_rows)
    }
    return [header] + new_rows, stats

# ----------------- Main -----------------

def main():
    if not os.path.exists(BASE_CELLS):
        print(f"Base cells.csv not found: {BASE_CELLS}")
        return 1
    header, data = read_cells(BASE_CELLS)
    z_open = find_third_layer_z(data)
    center = compute_center(data)
    print(f"Detected third layer z_open = {z_open}")
    print(f"Computed gland center (bounding box midpoint) = {center}")
    rows, stats = process(header, data, z_open, center)
    with open(OUT_SINGLE, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerows(rows)
    print("Single gland stats:", stats)
    print(f"Wrote: {OUT_SINGLE}")
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
