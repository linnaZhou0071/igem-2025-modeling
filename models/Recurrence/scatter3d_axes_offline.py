#!/usr/bin/env python3
"""
Improved fully-offline 3D scatter with axes / ticks / bounding box (no WebGL, no external libs).
Features:
  - Orthographic projection (default) for undistorted geometry.
  - Optionally perspective (--perspective) if you want depth.
  - Axis lines with tick marks + numeric labels (deg, aipc, vhh) and color legend.
  - Top percentile highlighting by larger white-ring points.
  - Stable equal-aspect scaling so the shape is not squashed.
  - Simple rotation by mouse drag; wheel zoom. Reset (key 'r'), orthographic/perspective toggle (key 'p').

Usage:
  python3 scatter3d_axes_offline.py --csv scored_runs.csv --top-percent 0.05 --output scatter3d_axes.html

If you want a static initial view only, just open the generated HTML (interactive is optional).
"""
import argparse, csv, math, json
from pathlib import Path

PLACEHOLDER_HTML = """<!DOCTYPE html><html><head><meta charset='utf-8'><title>Offline 3D Scatter (Axes)</title>
<style>
body { margin:0; font-family: sans-serif; background:#111; color:#e0e0e0; }
#panel { position:fixed; top:0; left:0; right:0; background:#1b1b1b; padding:6px 10px; font-size:13px; line-height:1.4; box-shadow:0 2px 4px rgba(0,0,0,0.4); }
#panel span { margin-right:14px; white-space:nowrap; }
#canvas { position:absolute; top:60px; left:0; width:100vw; height:calc(100vh - 60px); background:#000; cursor:grab; }
#tooltip { position:fixed; pointer-events:none; background:rgba(0,0,0,0.78); color:#fff; padding:6px 8px; font-size:12px; border:1px solid #666; border-radius:4px; display:none; white-space:pre; backdrop-filter:blur(2px); }
small { opacity:0.8; }
</style>
</head><body>
<div id='panel'>
  <span><b>3D 参数分布 (Axes)</b></span>
  <span>Points: __N__ (deg uniq: __U_DEG__, aipc uniq: __U_AIPC__, vhh uniq: __U_VHH__)</span>
  <span>Top%: __TOP__</span>
  <span>拖拽旋转 / 滚轮缩放 / r 重置 / p 投影切换 / e 模式切换(cube↔legacy) / g 稠密点 / d 调试</span>
  <span>模式: <span id='modeLabel'>Orthographic</span></span>
</div>
<canvas id='canvas'></canvas>
<div id='tooltip'></div>
<script>
const DATA = __DATA__;
const DENSE_DATA = __DENSE__; // optional interpolated grid (normalized 0..1 coords)
const TOP_IDS = new Set(__TOPIDS__);
// Extract arrays
let xs=DATA.map(d=>d.deg), ys=DATA.map(d=>d.aipc), zs=DATA.map(d=>d.vhh), scores=DATA.map(d=>d.Score);
// Unique values (for debug and verification)
let uniqDeg=[...new Set(xs)].sort((a,b)=>a-b);
let uniqAipc=[...new Set(ys)].sort((a,b)=>a-b);
let uniqVhh=[...new Set(zs)].sort((a,b)=>a-b);
let run_ids=DATA.map(d=>d.run_id), ranks=DATA.map(d=>d.rank);
function range(a){ return [Math.min(...a), Math.max(...a)]; }
const [xmin,xmax]=range(xs), [ymin,ymax]=range(ys), [zmin,zmax]=range(zs), [smin,smax]=range(scores);
// Cube normalization: map each axis linearly min->0, max->1
const spanX=(xmax-xmin)||1, spanY=(ymax-ymin)||1, spanZ=(zmax-zmin)||1;
// Raw (0..1) coordinates (origin at cube vertex = min values)
let RX=[],RY=[],RZ=[]; // raw normalized
for (let i=0;i<xs.length;i++){ RX.push((xs[i]-xmin)/spanX); RY.push((ys[i]-ymin)/spanY); RZ.push((zs[i]-zmin)/spanZ); }
// Internal centered coords for rotation: subtract 0.5
let X=[],Y=[],Z=[]; recomputeInternal();
function recomputeInternal(){
  X.length=0;Y.length=0;Z.length=0;
  for (let i=0;i<RX.length;i++){ X.push(RX[i]-0.5); Y.push(RY[i]-0.5); Z.push(RZ[i]-0.5); }
}
// Equal cube by construction; keep flag & key 'e' to toggle between corner-origin cube and legacy centered symmetric range if desired.
let cornerOrigin=true; // current requested mode
function toggleMode(){
  cornerOrigin=!cornerOrigin;
  if(cornerOrigin){ // rebuild from raw 0..1
    recomputeInternal();
  } else {
    // legacy: recenter but scale each axis to common max span (approx symmetrical) (mostly for comparison)
    const maxSpan=Math.max(spanX,spanY,spanZ);
    X.length=0;Y.length=0;Z.length=0;
    for (let i=0;i<xs.length;i++){
      X.push( (xs[i]- (xmin+xmax)/2)/maxSpan* (maxSpan/spanX) );
      Y.push( (ys[i]- (ymin+ymax)/2)/maxSpan* (maxSpan/spanY) );
      Z.push( (zs[i]- (zmin+zmax)/2)/maxSpan* (maxSpan/spanZ) );
    }
  }
}

// Viridis-like gradient
function viridis(t){ t=Math.max(0,Math.min(1,t)); const stops=[[68,1,84],[58,82,139],[33,145,140],[94,201,98],[170,220,50],[253,231,37]]; const idx=t*(stops.length-1); const i0=Math.floor(idx), i1=Math.min(stops.length-1,i0+1); const lt=idx-i0; const a=stops[i0],b=stops[i1]; const r=a[0]+(b[0]-a[0])*lt, g=a[1]+(b[1]-a[1])*lt, bl=a[2]+(b[2]-a[2])*lt; return `rgb(${r|0},${g|0},${bl|0})`; }
function norm(v,a,b){ return (v-a)/(b-a||1); }

const canvas=document.getElementById('canvas');
const ctx=canvas.getContext('2d');
let W,H; function resize(){ W=canvas.width=window.innerWidth; H=canvas.height=window.innerHeight-48; } window.addEventListener('resize',resize); resize();
let rotX=0.55, rotY=0.75; // nicer initial angle
let zoom=480; // slightly larger
let perspective=false; // start orthographic for geometric fidelity
let dragging=false,lastX=0,lastY=0;
canvas.addEventListener('mousedown',e=>{dragging=true;lastX=e.clientX;lastY=e.clientY;canvas.style.cursor='grabbing';});
window.addEventListener('mouseup',()=>{dragging=false;canvas.style.cursor='grab';});
window.addEventListener('mousemove',e=>{ if(dragging){ const dx=e.clientX-lastX, dy=e.clientY-lastY; rotY+=dx*0.005; rotX+=dy*0.005; lastX=e.clientX; lastY=e.clientY; requestRender(); }});
canvas.addEventListener('wheel',e=>{ e.preventDefault(); const f=Math.exp(-e.deltaY*0.001); zoom*=f; zoom=Math.max(40,Math.min(4000,zoom)); requestRender(); });
let debugOverlay=true;
let showDense=true; // toggle with 'g'
window.addEventListener('keydown',e=>{ if(e.key==='r'){rotX=0.55;rotY=0.75;zoom=480;requestRender();}
 else if(e.key==='p'){perspective=!perspective; document.getElementById('modeLabel').textContent=perspective?'Perspective':'Orthographic'; requestRender();}
 else if(e.key==='e'){toggleMode(); requestRender();}
 else if(e.key==='g'){showDense=!showDense; requestRender();}
 else if(e.key==='d'){debugOverlay=!debugOverlay; requestRender();}});

function project(x,y,z){ // rotate then project
  const cy=Math.cos(rotY), sy=Math.sin(rotY); const cx=Math.cos(rotX), sx=Math.sin(rotX);
  // rotate Y
  let xr=x*cy+z*sy; let zr=-x*sy+z*cy; // rotate X
  let yr=y*cx - zr*sx; let zr2=y*sx + zr*cx;
  if(perspective){ const d=1/(1+zr2*0.8); return {px:xr*zoom*d+W/2, py:-yr*zoom*d+H/2, depth:zr2}; }
  return {px:xr*zoom+W/2, py:-yr*zoom+H/2, depth:zr2};
}

let screenPts=[]; // cache
function render(){
  ctx.fillStyle='#000'; ctx.fillRect(0,0,W,H);
  // draw axes & bounding box first (behind points)
  drawBoxAndAxes();
  if(showDense && DENSE_DATA.length){
    // draw dense interpolated points (lighter / semi-transparent)
    const alphaBase=0.35; // transparency for dense layer
    for(const dp of DENSE_DATA){
      const pr=project(dp.x-0.5, dp.y-0.5, dp.z-0.5);
      const t=norm(dp.Score,smin,smax); const col=viridis(t);
      ctx.beginPath(); ctx.fillStyle=col.replace('rgb','rgba').replace(')',','+alphaBase+')');
      ctx.arc(pr.px,pr.py,2,0,Math.PI*2); ctx.fill();
    }
  }
  screenPts.length=0;
  for(let i=0;i<X.length;i++){
    const pr=project(X[i],Y[i],Z[i]); const c=viridis(norm(scores[i],smin,smax));
    screenPts.push({i,px:pr.px,py:pr.py,depth:pr.depth,color:c});
  }
  // painter sort (back to front) using depth ascending
  screenPts.sort((a,b)=>a.depth-b.depth);
  for(const p of screenPts){ const hi=TOP_IDS.has(run_ids[p.i]); ctx.beginPath(); ctx.fillStyle=p.color; ctx.arc(p.px,p.py, hi?5:3,0,Math.PI*2); ctx.fill(); if(hi){ ctx.strokeStyle='#fff'; ctx.lineWidth=1; ctx.stroke(); } }
  drawLegend();
}

function drawBoxAndAxes(){
  // Cube vertices in raw space (0..1) then convert to centered coords (subtract 0.5) for rotation
  const verts=[]; for (let dx of [0,1]) for (let dy of [0,1]) for (let dz of [0,1]) verts.push([dx-0.5,dy-0.5,dz-0.5]);
  // edges referencing binary-coded indices (dx + 2*dy + 4*dz) ordering
  const idx=(dx,dy,dz)=> dx + 2*dy + 4*dz;
  const edges=[
    [idx(0,0,0),idx(1,0,0)],[idx(1,0,0),idx(1,1,0)],[idx(1,1,0),idx(0,1,0)],[idx(0,1,0),idx(0,0,0)], // bottom
    [idx(0,0,1),idx(1,0,1)],[idx(1,0,1),idx(1,1,1)],[idx(1,1,1),idx(0,1,1)],[idx(0,1,1),idx(0,0,1)], // top
    [idx(0,0,0),idx(0,0,1)],[idx(1,0,0),idx(1,0,1)],[idx(1,1,0),idx(1,1,1)],[idx(0,1,0),idx(0,1,1)] // verticals
  ];
  ctx.strokeStyle='#2d2d2d'; ctx.lineWidth=1; ctx.beginPath();
  for(const [a,b] of edges){ const pa=project(...verts[a]); const pb=project(...verts[b]); ctx.moveTo(pa.px,pa.py); ctx.lineTo(pb.px,pb.py);} ctx.stroke();
  // axes (origin at cube corner raw (0,0,0) -> centered (-0.5,-0.5,-0.5)) to (1,0,0),(0,1,0),(0,0,1)
  const fullNames={deg:'vHH_to_cytokine_degradation_rate', aipc:'AIP_C_secretion_saturation', vhh:'vHH_secretion_saturation'};
  const axes=[ [ [-0.5,-0.5,-0.5],[0.65,-0.5,-0.5],'deg',[xmin,xmax] ], [ [-0.5,-0.5,-0.5],[-0.5,0.65,-0.5],'aipc',[ymin,ymax] ], [ [-0.5,-0.5,-0.5],[-0.5,-0.5,0.65],'vhh',[zmin,zmax] ] ];
  ctx.lineWidth=3; const axisColors={'deg':'#ff6b55','aipc':'#6bff55','vhh':'#55a8ff'};
  ctx.font='18px sans-serif'; ctx.textAlign='center'; ctx.textBaseline='middle';
  for(const [p0,p1,label,rangeVals] of axes){
    ctx.strokeStyle=axisColors[label]; ctx.beginPath(); const a=project(...p0), b=project(...p1); ctx.moveTo(a.px,a.py); ctx.lineTo(b.px,b.py); ctx.stroke();
    // arrow head
    const arrowSize=8; const dirX=(b.px-a.px), dirY=(b.py-a.py); const len=Math.hypot(dirX,dirY)||1; const ux=dirX/len, uy=dirY/len; ctx.beginPath(); ctx.moveTo(b.px,b.py); ctx.lineTo(b.px-ux*arrowSize-uy*arrowSize*0.6, b.py-uy*arrowSize+ux*arrowSize*0.6); ctx.lineTo(b.px-ux*arrowSize+uy*arrowSize*0.6, b.py-uy*arrowSize-ux*arrowSize*0.6); ctx.closePath(); ctx.fillStyle=axisColors[label]; ctx.fill();
    // label
    ctx.fillStyle=axisColors[label]; ctx.fillText(fullNames[label], b.px + 12, b.py);
    // ticks + grid planes
    const [vmin,vmax]=rangeVals; const ticks=10; for(let t=0;t<=ticks;t++){
      const frac=t/ticks; // raw tick in [0,1]
      let rawx= (label==='deg')? frac : 0;
      let rawy= (label==='aipc')? frac : 0;
      let rawz= (label==='vhh')? frac : 0;
      const vx = rawx - 0.5; const vy = rawy - 0.5; const vz = rawz - 0.5; // centered for projection
      const pt=project(vx,vy,vz);
      // tick mark
      ctx.beginPath(); ctx.fillStyle='#ddd'; ctx.arc(pt.px,pt.py,3.5,0,Math.PI*2); ctx.fill();
      // text slightly offset
      const rawVal=(vmin + (vmax-vmin)*frac);
      ctx.fillStyle='#bbb'; ctx.font='14px sans-serif'; ctx.fillText(rawVal.toFixed(2), pt.px, pt.py - 16);
    }
  }
  if(debugOverlay){ drawDebug(); }
}

function drawLegend(){ const w=230,h=18,x=24,y=24; ctx.save(); ctx.globalAlpha=0.95; for(let i=0;i<w;i++){ const t=i/(w-1); ctx.fillStyle=viridis(t); ctx.fillRect(x+i,y,1,h);} ctx.restore(); ctx.strokeStyle='#aaa'; ctx.strokeRect(x,y,w,h); ctx.fillStyle='#eee'; ctx.font='14px sans-serif'; ctx.fillText(smin.toFixed(3), x, y+h+18); ctx.fillText(smax.toFixed(3), x+w-50, y+h+18); ctx.fillText('Score', x+w/2-30, y-6); }

function drawDebug(){
  const pad=18; let lines=[
    'deg range: '+xmin.toFixed(3)+' – '+xmax.toFixed(3)+' (uniq '+uniqDeg.length+')',
    'aipc range: '+ymin.toFixed(3)+' – '+ymax.toFixed(3)+' (uniq '+uniqAipc.length+')',
    'vhh range: '+zmin.toFixed(3)+' – '+zmax.toFixed(3)+' (uniq '+uniqVhh.length+')',
    'mode=cubeOrigin perspective='+perspective,
  ];
  ctx.save(); ctx.font='13px monospace'; ctx.textAlign='left'; ctx.textBaseline='top';
  const boxW=360; const boxH=lines.length*18 + 10; ctx.globalAlpha=0.55; ctx.fillStyle='#000'; ctx.fillRect(W-boxW-pad, pad, boxW, boxH); ctx.globalAlpha=1; ctx.strokeStyle='#555'; ctx.strokeRect(W-boxW-pad, pad, boxW, boxH);
  let y0=pad+6; ctx.fillStyle='#8fd1ff'; for(const ln of lines){ ctx.fillText(ln, W-boxW-pad+10, y0); y0+=18; }
  ctx.restore();
}

let needRender=true; function requestRender(){ needRender=true; }
function loop(){ if(needRender){ render(); needRender=false; } requestAnimationFrame(loop);} loop(); requestRender();

// Tooltip
const tooltip=document.getElementById('tooltip');
canvas.addEventListener('mousemove',e=>{ const mx=e.clientX, my=e.clientY-48; let best=null, bestd=100; for(const p of screenPts){ const dx=p.px-mx, dy=p.py-my; const d=dx*dx+dy*dy; if(d<bestd){bestd=d;best=p;} } if(best && bestd<100){ const i=best.i; tooltip.style.display='block'; tooltip.style.left=(e.clientX+12)+'px'; tooltip.style.top=(e.clientY+12)+'px'; tooltip.textContent=`run ${run_ids[i]}\nrank ${ranks[i]}\nScore=${scores[i].toFixed(4)}\ndeg=${xs[i]}\naipc=${ys[i]}\nvhh=${zs[i]}`; } else { tooltip.style.display='none'; } });
</script>
</body></html>"""

def read_rows(path: Path):
    rows=[]
    with path.open('r', newline='') as f:
        r=csv.DictReader(f)
        for row in r:
            try:
                rows.append({
                    'Score': float(row['Score']),
                    'run_id': row.get('run_id',''),
                    'rank': int(row.get('rank','0') or 0),
                    'deg': float(row['degradation_rate']),
                    'aipc': float(row['aipc_saturation']),
                    'vhh': float(row['vhh_saturation'])
                })
            except Exception:
                continue
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--csv', type=Path, default=Path('scored_runs.csv'))
    ap.add_argument('--top-percent', type=float, default=0.05)
    ap.add_argument('--output', type=Path, default=Path('scatter3d_axes.html'))
    ap.add_argument('--dense', type=int, default=0, help='Per-axis dense interpolation resolution (e.g. 25 -> 25^3 points). 0 disables.')
    args=ap.parse_args()

    rows=read_rows(args.csv)
    if not rows:
        print('[ERROR] no rows read')
        return
    has_rank=all('rank' in r for r in rows)
    if has_rank:
        rows_sorted=sorted(rows, key=lambda r: r['rank'])
    else:
        rows_sorted=sorted(rows, key=lambda r: r['Score'], reverse=True)
    top_n=max(1,int(len(rows)*args.top_percent))
    top_ids=[r['run_id'] for r in rows_sorted[:top_n]]
    # Optional dense interpolation grid (trilinear on structured 10x10x10 lattice)
    dense_list = []
    if args.dense and args.dense > 1:
        if args.dense > 40:
            print('[warn] --dense capped at 40 for performance')
            args.dense = 40
        deg_vals = sorted({r['deg'] for r in rows})
        aipc_vals = sorted({r['aipc'] for r in rows})
        vhh_vals = sorted({r['vhh'] for r in rows})
        idx_deg = {v: i for i, v in enumerate(deg_vals)}
        idx_aipc = {v: i for i, v in enumerate(aipc_vals)}
        idx_vhh = {v: i for i, v in enumerate(vhh_vals)}
        cube = [[[0.0 for _ in vhh_vals] for _ in aipc_vals] for _ in deg_vals]
        for r in rows:
            cube[idx_deg[r['deg']]][idx_aipc[r['aipc']]][idx_vhh[r['vhh']]] = r['Score']
        n = args.dense
        for ii in range(n):
            u = ii / (n - 1) * (len(deg_vals) - 1)
            i0 = int(math.floor(u)); i1 = min(i0 + 1, len(deg_vals) - 1); du = u - i0
            if i0 == i1: du = 0.0
            for jj in range(n):
                v = jj / (n - 1) * (len(aipc_vals) - 1)
                j0 = int(math.floor(v)); j1 = min(j0 + 1, len(aipc_vals) - 1); dv = v - j0
                if j0 == j1: dv = 0.0
                for kk in range(n):
                    w = kk / (n - 1) * (len(vhh_vals) - 1)
                    k0 = int(math.floor(w)); k1 = min(k0 + 1, len(vhh_vals) - 1); dw = w - k0
                    if k0 == k1: dw = 0.0
                    c000 = cube[i0][j0][k0]; c100 = cube[i1][j0][k0]; c010 = cube[i0][j1][k0]; c110 = cube[i1][j1][k0]
                    c001 = cube[i0][j0][k1]; c101 = cube[i1][j0][k1]; c011 = cube[i0][j1][k1]; c111 = cube[i1][j1][k1]
                    c00 = c000 * (1 - du) + c100 * du; c01 = c001 * (1 - du) + c101 * du
                    c10 = c010 * (1 - du) + c110 * du; c11 = c011 * (1 - du) + c111 * du
                    c0 = c00 * (1 - dv) + c10 * dv; c1 = c01 * (1 - dv) + c11 * dv
                    c = c0 * (1 - dw) + c1 * dw
                    dense_list.append({
                        'x': (deg_vals[0] + (deg_vals[-1] - deg_vals[0]) * (ii / (n - 1)) - deg_vals[0]) / (deg_vals[-1] - deg_vals[0] or 1),
                        'y': (aipc_vals[0] + (aipc_vals[-1] - aipc_vals[0]) * (jj / (n - 1)) - aipc_vals[0]) / (aipc_vals[-1] - aipc_vals[0] or 1),
                        'z': (vhh_vals[0] + (vhh_vals[-1] - vhh_vals[0]) * (kk / (n - 1)) - vhh_vals[0]) / (vhh_vals[-1] - vhh_vals[0] or 1),
                        'Score': float(c)
                    })
        print(f'[info] dense grid: {n}^3 = {n**3} points')
    html=(PLACEHOLDER_HTML
          .replace('__N__', str(len(rows)))
          .replace('__TOP__', str(args.top_percent))
          .replace('__DATA__', json.dumps(rows))
          .replace('__TOPIDS__', json.dumps(top_ids))
          .replace('__DENSE__', json.dumps(dense_list)))
    args.output.write_text(html)
    print('[done] wrote', args.output)

if __name__=='__main__':
    main()
