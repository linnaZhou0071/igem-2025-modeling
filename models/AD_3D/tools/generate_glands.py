#!/usr/bin/env python3
"""
Generate modified cells.csv with two sebaceous gland geometries (variants A and D).

Assumptions / Notes:
 - Base file: user_projects/AD_3D/config/cells.csv
 - Only keratinocyte rows are modified / removed / added. All other cell types are preserved.
 - Two gland centers are fixed (deterministic): (-120, 80) and (130, -60)
 - Variant A: radius expands linearly from 50 -> 150 over 150 um depth. Bottom closed by adding a bottom disc (same radius) plus one inner ring (approximating a shallow cap) without extending depth beyond 150 um.
 - Variant D: neck (first 40 um) fixed radius 50; then smooth (power 1.2) expansion to 150 at 150 um depth.
 - Mild bending (amplitude 10 um): gland 1 bends in +x via sin(pi*d/total_depth); gland 2 bends in +y.
 - Third layer (z_open) determined by rounding keratinocyte z to nearest integer, sorting unique z descending, and taking the 3rd element.
 - Any keratinocyte within the gland volume (radial distance < current radius - shell_margin) is removed (H1).
 - Output files: cells_with_glands_A.csv, cells_with_glands_D.csv placed beside original file.
"""

from __future__ import annotations
import csv, math, os, sys
from collections import defaultdict

BASE_CELLS = os.path.join(os.path.dirname(__file__), "..", "config", "cells.csv")
OUT_A = os.path.join(os.path.dirname(BASE_CELLS), "cells_with_glands_A.csv")
OUT_D = os.path.join(os.path.dirname(BASE_CELLS), "cells_with_glands_D.csv")

GlandCenters = [(-120.0, 80.0), (130.0, -60.0)]  # (cx, cy)
TopRadius = 50.0
BottomRadius = 150.0
TotalDepth = 150.0
LayerStep = 10.0  # axial discretization (um)
ShellMargin = 4.0 # how far inside the computed radius we purge existing keratinocytes
CircumSpacing = 15.0  # target arc length spacing (variant F2)
BendAmp = 10.0

def read_cells(path):
    with open(path) as f:
        reader = csv.reader(f)
        rows = list(reader)
    header = rows[0]
    data = rows[1:]
    return header, data

def parse_row(row):
    # x,y,z,type,... -> floats for coords
    x = float(row[0]); y = float(row[1]); z = float(row[2]); ctype = row[3]
    return x,y,z,ctype

def find_third_layer_z(data):
    z_vals = []
    for r in data:
        try:
            x,y,z,ctype = parse_row(r)
        except Exception:
            continue
        if ctype.strip() != 'keratinocyte':
            continue
        z_vals.append(round(z))  # B2 rounding
    uniq = sorted(set(z_vals), reverse=True)
    if len(uniq) < 3:
        raise RuntimeError(f"Not enough distinct keratinocyte z-layers (found {uniq})")
    return float(uniq[2])

def centerline_offset(g_index, depth):
    # bending E2: gland1 bends in x, gland2 in y
    phase = math.pi * depth / TotalDepth
    dx = BendAmp * math.sin(phase) if g_index == 0 else 0.0
    dy = BendAmp * math.sin(phase) if g_index == 1 else 0.0
    return dx, dy

def radius_A(depth):
    # linear expansion 50->150 over 0..150
    return TopRadius + (BottomRadius - TopRadius) * (depth / TotalDepth)

def radius_D(depth):
    # neck first 40 um fixed 50; then smooth power ramp to 150
    if depth <= 40.0:
        return TopRadius
    frac = (depth - 40.0) / (TotalDepth - 40.0)
    return TopRadius + (BottomRadius - TopRadius) * (frac ** 1.2)

def build_shell(z_open, variant):
    radii_func = radius_A if variant == 'A' else radius_D
    shell_points = []  # (x,y,z)
    depths = [d for d in [i*LayerStep for i in range(int(TotalDepth/LayerStep)+1)] if d <= TotalDepth]
    for gi,(cx,cy) in enumerate(GlandCenters):
        for d in depths:
            r = radii_func(d)
            if r <= 0: continue
            dx_c, dy_c = centerline_offset(gi, d)
            cx_d = cx + dx_c
            cy_d = cy + dy_c
            # circumferential discretization
            n_theta = max(6, round(2*math.pi*r / CircumSpacing))
            for k in range(n_theta):
                theta = 2*math.pi*k / n_theta
                x = cx_d + r*math.cos(theta)
                y = cy_d + r*math.sin(theta)
                z = z_open - d
                shell_points.append((x,y,z))
        # bottom closure (simple: disc + inner ring) at d=TotalDepth
        rB = radii_func(TotalDepth)
        dx_c, dy_c = centerline_offset(gi, TotalDepth)
        cx_b = cx + dx_c; cy_b = cy + dy_c; z_b = z_open - TotalDepth
        n_theta = max(8, round(2*math.pi*rB / CircumSpacing))
        # inner ring at 0.6*radius to approximate curved closure
        for ring_frac in (0.6, 0.3):
            rr = rB * ring_frac
            n_ring = max(6, round(2*math.pi*rr / CircumSpacing))
            for k in range(n_ring):
                theta = 2*math.pi*k / n_ring
                x = cx_b + rr*math.cos(theta)
                y = cy_b + rr*math.sin(theta)
                shell_points.append((x,y,z_b))
        # center point
        shell_points.append((cx_b, cy_b, z_b))
    return shell_points

def inside_gland_volume(x,y,z,z_open,variant):
    depth = z_open - z
    if depth < 0 or depth > TotalDepth: return False
    rf = radius_A if variant=='A' else radius_D
    for gi,(cx,cy) in enumerate(GlandCenters):
        dx_c, dy_c = centerline_offset(gi, depth)
        cx_d = cx + dx_c; cy_d = cy + dy_c
        r_allow = rf(depth) - ShellMargin
        if r_allow <= 0: continue
        if (x-cx_d)**2 + (y-cy_d)**2 < r_allow**2:
            return True
    return False

def cylinder_cut(x,y):
    for (cx,cy) in GlandCenters:
        if (x-cx)**2 + (y-cy)**2 < TopRadius**2:
            return True
    return False

def process_variant(header, data, variant, z_open):
    new_rows = []
    cut_removed = 0
    vol_removed = 0
    for row in data:
        try:
            x,y,z,ctype = parse_row(row)
        except Exception:
            continue
        if ctype == 'keratinocyte':
            # initial cylinder cut
            if cylinder_cut(x,y):
                cut_removed += 1
                continue
            # gland volume purge
            if inside_gland_volume(x,y,z,z_open,variant):
                vol_removed += 1
                continue
        new_rows.append(row)
    # add shell cells
    shell = build_shell(z_open, variant)
    for (x,y,z) in shell:
        new_rows.append([f"{x:.3f}", f"{y:.3f}", f"{z:.3f}", 'keratinocyte', '', '', '', ''])
    stats = {
        'cut_removed': cut_removed,
        'vol_removed': vol_removed,
        'shell_added': len(shell),
        'final_total': len(new_rows)
    }
    return [header] + new_rows, stats

def main():
    if not os.path.exists(BASE_CELLS):
        print(f"Base cells.csv not found: {BASE_CELLS}")
        sys.exit(1)
    header, data = read_cells(BASE_CELLS)
    z_open = find_third_layer_z(data)
    print(f"Detected third layer z_open = {z_open}")
    # Variant A
    rowsA, statsA = process_variant(header, data, 'A', z_open)
    with open(OUT_A, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(rowsA)
    # Variant D
    rowsD, statsD = process_variant(header, data, 'D', z_open)
    with open(OUT_D, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(rowsD)
    print("Variant A stats:", statsA)
    print("Variant D stats:", statsD)
    print(f"Wrote: {OUT_A}\n       {OUT_D}")

if __name__ == '__main__':
    main()
