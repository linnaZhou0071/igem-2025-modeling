#!/usr/bin/env python3
import csv,math,random,os,sys

CSV_PATH = os.path.join(os.path.dirname(__file__),"cells.csv")
BACKUP_PATH = os.path.join(os.path.dirname(__file__),"cells_pre_new_domes_backup.csv")
NEW_TYPE = "engr.S.epidermis"
NUM_DOMES = 10
CELLS_PER_DOME = 50
BASE_Z = 260.0             # keratinocyte surface plane
SEED = 24680
PACKING_FRACTION = 0.60
ENG_DIAMETER = 7.0
ENG_RADIUS = ENG_DIAMETER/2.0
ENG_CELL_VOLUME = (4.0/3.0)*math.pi*(ENG_RADIUS**3)
RADIUS_SCALE = 1.05        # slight looseness vs densest theoretical
def colony_radius(n):
    target_vol = n*ENG_CELL_VOLUME / PACKING_FRACTION
    # hemisphere volume = 2/3 pi R^3
    R3 = target_vol / ((2.0/3.0)*math.pi)
    return (R3**(1.0/3.0))*RADIUS_SCALE
MIN_CENTER_DIST = 2.5 * ENG_DIAMETER * (CELLS_PER_DOME/50.0)**(1/3)
# Sebaceous gland void (exclude centers & points inside). You can tune these to match AD_3D gland if needed.
# Here assume gland roughly centered at (0,0) with outer radius ~80 and inner open column ~40 at top.
GLAND_CX = 0.0
GLAND_CY = 0.0
GLAND_EXCLUDE_RADIUS = 50.0  # do not place dome centers inside this radius (approx opening)

random.seed(SEED)

if not os.path.isfile(CSV_PATH):
    print("Missing cells.csv")
    sys.exit(1)

with open(CSV_PATH,'r',newline='') as f:
    rows = list(csv.reader(f))
header = rows[0]
body = rows[1:]

# Remove old engr.S.epidermis if any
filtered=[]; removed=0
for r in body:
    if len(r)>=4 and r[3].strip()==NEW_TYPE:
        removed+=1
    else:
        filtered.append(r)
body=filtered

# Determine (x,y) bounds from keratinocytes at z near BASE_Z- some tolerance
ker_points=[]
for r in body:
    if len(r)>=4 and r[3].strip()=="keratinocyte":
        try:
            x=float(r[0]); y=float(r[1]); z=float(r[2])
        except: continue
        if abs(z-BASE_Z) < 10: # top layer band
            ker_points.append((x,y))

if not ker_points:
    print("No keratinocyte points near z=260; cannot anchor domes.")
    sys.exit(1)
xs=[p[0] for p in ker_points]; ys=[p[1] for p in ker_points]
XMIN,XMAX=min(xs),max(xs); YMIN,YMAX=min(ys),max(ys)

# Place centers avoiding gland opening and each other
centers=[]; attempts=0
while len(centers)<NUM_DOMES and attempts<8000:
    attempts+=1
    cx=random.uniform(XMIN,XMAX)
    cy=random.uniform(YMIN,YMAX)
    # exclude gland opening
    if (cx-GLAND_CX)**2 + (cy-GLAND_CY)**2 < GLAND_EXCLUDE_RADIUS**2:
        continue
    ok=True
    for (ox,oy) in centers:
        if math.hypot(cx-ox,cy-oy) < MIN_CENTER_DIST:
            ok=False; break
    if ok:
        centers.append((cx,cy))
if len(centers)<NUM_DOMES:
    print(f"WARNING: only placed {len(centers)} domes")

def generate_hemisphere_points(n, R, attempts_factor=600):
    pts=[]
    min_sep = ENG_DIAMETER*0.85
    min_sep2 = min_sep*min_sep
    max_attempts = n * attempts_factor
    attempts=0
    while len(pts) < n and attempts < max_attempts:
        attempts+=1
        px = random.uniform(-R,R)
        py = random.uniform(-R,R)
        pz = random.uniform(0,R)
        if px*px + py*py + pz*pz > R*R: continue
        good=True
        for (ax,ay,az) in pts:
            dx=px-ax; dy=py-ay; dz=pz-az
            if dx*dx+dy*dy+dz*dz < min_sep2: good=False; break
        if not good: continue
        pts.append((px,py,pz))
    return pts, attempts

new_rows=[]
colony_info=[]
for ci,(cx,cy) in enumerate(centers):
    R = colony_radius(CELLS_PER_DOME)
    pts, atts = generate_hemisphere_points(CELLS_PER_DOME, R)
    if len(pts) < CELLS_PER_DOME:
        print(f"WARNING colony {ci}: only {len(pts)}/{CELLS_PER_DOME} placed")
    for (px,py,pz) in pts:
        new_rows.append([f"{cx+px:.3f}", f"{cy+py:.3f}", f"{BASE_Z+pz:.3f}", NEW_TYPE, "","","",""])
    colony_info.append((ci,R,len(pts),atts))

# Backup original (without previous engr)
if not os.path.exists(BACKUP_PATH):
    with open(BACKUP_PATH,'w',newline='') as bf:
        w=csv.writer(bf); w.writerow(header); w.writerows(body)

# Write updated CSV
with open(CSV_PATH,'w',newline='') as f:
    w=csv.writer(f); w.writerow(header); w.writerows(body); w.writerows(new_rows)

print(f"Removed old {NEW_TYPE}: {removed}")
print(f"Placed {len(centers)} colony centers (hemispherical random packing)")
print(f"Added cells: {len(new_rows)} (expected {NUM_DOMES*CELLS_PER_DOME})")
for (ci,R,placed,atts) in colony_info:
    print(f"  Colony {ci}: R={R:.2f} cells={placed} attempts={atts}")
print("Centers:")
for i,(cx,cy) in enumerate(centers):
    print(f"  {i}: ({cx:.2f},{cy:.2f})")
print("DONE")
