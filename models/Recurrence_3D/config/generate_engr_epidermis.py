#!/usr/bin/env python3
import csv,random,math,os,sys

CSV_PATH = os.path.join(os.path.dirname(__file__),"cells.csv")
BACKUP_PATH = os.path.join(os.path.dirname(__file__),"cells_pre_engr_backup.csv")
NEW_TYPE = "engr.S.epidermis"
NUM_COLONIES = 10
CELLS_PER_COLONY = 50
COLONY_RADIUS = 30.0      # microns (approx lateral radius)
FIX_TOP_Z = 260.0         # desired apex z
BASE_OFFSET = 0.5         # small lift above the keratinocyte top plane to keep colony edges from coinciding with it
MIN_CENTER_SEPARATION = 0.6 * 2 * COLONY_RADIUS  # allow some overlap but avoid heavy collision
RANDOM_SEED = 12345

if not os.path.isfile(CSV_PATH):
    print(f"ERROR: {CSV_PATH} not found")
    sys.exit(1)

# Read existing
with open(CSV_PATH,'r',newline='') as f:
    rows = list(csv.reader(f))

if not rows:
    print("Empty CSV")
    sys.exit(1)
header = rows[0]
body = rows[1:]

# Remove previously generated engr.S.epidermis (idempotent regeneration)
filtered_body = []
removed_prev = 0
for r in body:
    if len(r) >= 4 and r[3].strip() == NEW_TYPE:
        removed_prev += 1
    else:
        filtered_body.append(r)
body = filtered_body

# Basic header validation
expected_header = ['x','y','z','type','volume','cycle entry','custom:GFP','custom:sample']
if header[:len(expected_header)] != expected_header:
    print("WARNING: CSV header unexpected, proceeding anyway:", header)

# Extract keratinocyte stats
ker_x=[]; ker_y=[]; ker_z=[]
for r in body:
    if len(r) < 4: continue
    if r[3].strip()=="keratinocyte":
        try:
            ker_x.append(float(r[0])); ker_y.append(float(r[1])); ker_z.append(float(r[2]))
        except ValueError:
            pass

if not ker_z:
    print("ERROR: No keratinocyte rows found to anchor colonies.")
    sys.exit(1)

x_min,x_max = min(ker_x), max(ker_x)
y_min,y_max = min(ker_y), max(ker_y)
z_base_plane = max(ker_z)  # detected top keratinocyte layer

# Dynamic dome height so that apex == FIX_TOP_Z (clamped >= minimal)
raw_height = FIX_TOP_Z - (z_base_plane + BASE_OFFSET)
if raw_height < 0.5:
    raw_height = 0.5  # ensure some positive thickness
DOME_HEIGHT = raw_height

# Place colony centers randomly in bounding box with spacing constraint
random.seed(RANDOM_SEED)
centers=[]
tries=0
while len(centers) < NUM_COLONIES and tries < 5000:
    tries += 1
    cx = random.uniform(x_min, x_max)
    cy = random.uniform(y_min, y_max)
    ok=True
    for (ox,oy) in centers:
        if math.hypot(cx-ox, cy-oy) < MIN_CENTER_SEPARATION:
            ok=False; break
    if ok:
        centers.append((cx,cy))

if len(centers) < NUM_COLONIES:
    print(f"WARNING: Only placed {len(centers)} colony centers (requested {NUM_COLONIES}).")

# Deterministic point generation inside each colony
new_rows=[]
for idx,(cx,cy) in enumerate(centers):
    for k in range(CELLS_PER_COLONY):
        # area-uniform radius
        u = ( (k+1)*(idx+3)*0.318309886 ) % 1.0   # frac using 1/pi
        v = ( (k+7)*(idx+11)*0.141421356 ) % 1.0  # frac using (sqrt(2)-1) approx
        r = COLONY_RADIUS * math.sqrt(u)
        theta = 2*math.pi * v
        x = cx + r*math.cos(theta)
        y = cy + r*math.sin(theta)
        # Parabolic dome: height = H*(1 - (r/R)^2)
        z = z_base_plane + 0.5 + DOME_HEIGHT * (1.0 - (r/COLONY_RADIUS)**2)
        # Format similar to existing rows: volume and custom fields blank
        new_rows.append([f"{x:.3f}", f"{y:.3f}", f"{z:.3f}", NEW_TYPE, "", "", "", ""])

print(f"Removed previous {NEW_TYPE} cells: {removed_prev}")
print(f"Base keratinocyte detected top z = {z_base_plane:.3f}")
print(f"Target apex z = {FIX_TOP_Z:.3f}; computed dome height = {DOME_HEIGHT:.3f}")
print(f"Generated {len(new_rows)} new {NEW_TYPE} cells across {len(centers)} colonies.")

# Backup original
if not os.path.exists(BACKUP_PATH):
    with open(BACKUP_PATH,'w',newline='') as bf:
        w = csv.writer(bf)
        w.writerow(header)
        w.writerows(body)
    print(f"Backup written to {BACKUP_PATH}")
else:
    print("Backup already exists, not overwriting:", BACKUP_PATH)

# Append and write back
with open(CSV_PATH,'w',newline='') as f:
    w = csv.writer(f)
    w.writerow(header)
    w.writerows(body)
    w.writerows(new_rows)
print(f"Updated file with total rows = {1+len(body)+len(new_rows)}")

# Simple verification summary
print("Verification summary:")
print("  Colony centers:")
for i,(cx,cy) in enumerate(centers):
    print(f"    {i}: ({cx:.1f},{cy:.1f})")
print(f"  Each colony target cells: {CELLS_PER_COLONY}")
print(f"  Distinct colonies placed: {len(centers)}")
