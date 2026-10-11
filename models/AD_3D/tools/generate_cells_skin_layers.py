#!/usr/bin/env python3
"""
Generate layered 3D initial conditions for AD_3D project:
 - 3 SA (S.aurenus) microcolonies (hemispherical) on top surface (epidermis surface at z = 260)
 - 3 layers of keratinocytes (epidermis) thickness 60 um (each layer 20 um) occupying z in [200,260)
 - Dermis thickness 200 um (z in [0,200)) with 3 T cell spherical clusters (50 cells each)
 - Diameters: SA 7 um, T cell 7 um, keratinocyte 20 um
 - Sebocyte ignored for now.
Result CSV columns: x,y,z,type,volume,cycle entry,custom:GFP,custom:sample
Volumes computed from diameter assuming sphere.
"""
import math, random, csv, os
random.seed(42)

# Parameters
AREA_SIZE = 500.0  # square side length (x,y span)
X_MIN = -AREA_SIZE/2
Y_MIN = -AREA_SIZE/2
X_MAX = AREA_SIZE/2
Y_MAX = AREA_SIZE/2
DERMIS_THICKNESS = 200.0
EPIDERMIS_THICKNESS = 60.0
TOP_SURFACE_Z = DERMIS_THICKNESS + EPIDERMIS_THICKNESS  # 260

# Cell type specs (diameter microns)
DIAM = {
    'S.aurenus': 7.0,
    'keratinocyte': 20.0,
    'T_cell': 7.0,
}

def volume_from_diameter(d):
    r = d/2.0
    return (4.0/3.0)*math.pi*(r**3)

VOLUME = {k: volume_from_diameter(v) for k,v in DIAM.items()}

# Microcolony specifications: (count, center_x, center_y)
# Radii will be inferred from count and cell volume (approx hemispherical packing at 60% packing efficiency)
MICROCOLONIES = [
    (20, -120.0, -120.0),
    (30, 130.0, 50.0),
    (50, -50.0, 140.0),
]
PACKING_FRACTION = 0.60  # heuristic effective fraction of volume filled by cells

def hemisphere_radius(n_cells, cell_type):
    cell_vol = VOLUME[cell_type]
    target_volume = n_cells * cell_vol / PACKING_FRACTION
    # Hemisphere volume = (2/3) pi r^3
    r_cubed = target_volume / ((2.0/3.0)*math.pi)
    return r_cubed ** (1.0/3.0)

# T cell cluster specs: (count, center_x, center_y, center_z)
T_CELL_CLUSTERS = [
    (50, -150.0, -150.0, 60.0),
    (50, 150.0, 0.0, 110.0),
    (50, 0.0, 150.0, 150.0),
]

def sphere_radius(n_cells, cell_type):
    cell_vol = VOLUME[cell_type]
    target_volume = n_cells * cell_vol / PACKING_FRACTION
    # Sphere volume = 4/3 pi r^3
    r_cubed = target_volume / ((4.0/3.0)*math.pi)
    return r_cubed ** (1.0/3.0)

# Keratinocyte grid layers
KER_LAYER_COUNT = 3
KER_LAYER_THICKNESS = EPIDERMIS_THICKNESS / KER_LAYER_COUNT  # 20
KER_Z_LEVELS = [DERMIS_THICKNESS + (i+0.5)*KER_LAYER_THICKNESS for i in range(KER_LAYER_COUNT)]  # midpoints
KER_SPACING = DIAM['keratinocyte']  # ~20 um spacing
JITTER = 2.0

cells = []  # list of dicts

# Generate keratinocyte layers
x_positions = []
y_positions = []
cur = X_MIN + DIAM['keratinocyte']/2.0
while cur <= X_MAX - DIAM['keratinocyte']/2.0 + 1e-6:
    x_positions.append(cur)
    cur += KER_SPACING
cur = Y_MIN + DIAM['keratinocyte']/2.0
while cur <= Y_MAX - DIAM['keratinocyte']/2.0 + 1e-6:
    y_positions.append(cur)
    cur += KER_SPACING

for zc in KER_Z_LEVELS:
    for x in x_positions:
        for y in y_positions:
            cells.append({
                'x': x + random.uniform(-JITTER, JITTER),
                'y': y + random.uniform(-JITTER, JITTER),
                'z': zc + random.uniform(-JITTER, JITTER),
                'type': 'keratinocyte',
                'volume': VOLUME['keratinocyte'],
            })

# Generate T cell clusters
for n, cx, cy, cz in T_CELL_CLUSTERS:
    r = sphere_radius(n, 'T_cell')
    count = 0
    attempts = 0
    while count < n and attempts < n*500:
        # random point in cube then reject outside sphere
        px = random.uniform(-r, r)
        py = random.uniform(-r, r)
        pz = random.uniform(-r, r)
        if px*px + py*py + pz*pz <= r*r:
            cells.append({
                'x': cx + px,
                'y': cy + py,
                'z': cz + pz,
                'type': 'T_cell',
                'volume': VOLUME['T_cell'],
            })
            count += 1
        attempts += 1

# Generate SA microcolonies (hemispheres above surface)
for n, cx, cy in MICROCOLONIES:
    r = hemisphere_radius(n, 'S.aurenus')
    count = 0
    attempts = 0
    while count < n and attempts < n*800:
        # sample in full sphere then keep z>=0 (hemisphere)
        px = random.uniform(-r, r)
        py = random.uniform(-r, r)
        pz = random.uniform(0, r)  # directly sample positive z to improve efficiency
        if px*px + py*py + pz*pz <= r*r:
            cells.append({
                'x': cx + px,
                'y': cy + py,
                'z': TOP_SURFACE_Z + pz,  # sits above surface
                'type': 'S.aurenus',
                'volume': VOLUME['S.aurenus'],
            })
            count += 1
        attempts += 1

# Write CSV
out_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'cells.csv')
# Normalize path
out_path = os.path.normpath(out_path)

with open(out_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['x','y','z','type','volume','cycle entry','custom:GFP','custom:sample'])
    for c in cells:
        w.writerow([f"{c['x']:.3f}", f"{c['y']:.3f}", f"{c['z']:.3f}", c['type'], f"{c['volume']:.3f}", '', '0', '0'])

print(f"Wrote {len(cells)} cells to {out_path}")
