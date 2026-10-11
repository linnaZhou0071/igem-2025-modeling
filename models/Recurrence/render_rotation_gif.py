#!/usr/bin/env python3
"""Render a rotating parameter cube (GIF + optional MP4) with enlarged labels.

Features:
  * Uses scored_runs.csv (columns: degradation_rate, aipc_saturation, vhh_saturation, Score, rank, run_id)
  * Normalizes the 3 parameters to a unit cube and rotates around selected axis.
  * Dense trilinear interpolation background (--dense N) for smooth color field.
  * Highlight top fraction (--highlight-top) with white ring.
  * 3x enlarged fonts / ticks / legend for presentation clarity.
  * Adjustable rotation coverage via --sweep-angle * --angle-scale (full 360° default).
  * Optional MP4 export using ffmpeg (--export-mp4).

Example:
  python3 render_rotation_gif.py --csv scored_runs.csv --output rotation_fast.gif \
      --frames 120 --fps 24 --axis y --dense 25 --highlight-top 0.05 \
      --angle-scale 1.0 --export-mp4
"""
from __future__ import annotations
import csv, math, argparse, sys, shutil, subprocess
from pathlib import Path
from typing import List, Dict
from PIL import Image, ImageDraw, ImageFont

def read_rows(path: Path) -> List[Dict]:
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

def viridis(t: float):
    t=max(0.0,min(1.0,t))
    stops=[(68,1,84),(58,82,139),(33,145,140),(94,201,98),(170,220,50),(253,231,37)]
    x=t*(len(stops)-1); i0=int(math.floor(x)); i1=min(len(stops)-1,i0+1); lt=x-i0
    a=stops[i0]; b=stops[i1]
    r=int(a[0]+(b[0]-a[0])*lt); g=int(a[1]+(b[1]-a[1])*lt); bl=int(a[2]+(b[2]-a[2])*lt)
    return (r,g,bl)

def trilinear_dense(rows: List[Dict], n: int):
    if n<=1: return []
    if n>40: print('[warn] dense >40 capped to 40'); n=40
    deg_vals=sorted({r['deg'] for r in rows}); aipc_vals=sorted({r['aipc'] for r in rows}); vhh_vals=sorted({r['vhh'] for r in rows})
    idx_deg={v:i for i,v in enumerate(deg_vals)}; idx_aipc={v:i for i,v in enumerate(aipc_vals)}; idx_vhh={v:i for i,v in enumerate(vhh_vals)}
    cube=[ [ [0.0 for _ in vhh_vals] for _ in aipc_vals] for _ in deg_vals ]
    for r in rows: cube[idx_deg[r['deg']]][idx_aipc[r['aipc']]][idx_vhh[r['vhh']]]=r['Score']
    dens=[]
    for ii in range(n):
        u=ii/(n-1)*(len(deg_vals)-1); i0=int(u); i1=min(i0+1,len(deg_vals)-1); du=0 if i0==i1 else u-i0
        for jj in range(n):
            v=jj/(n-1)*(len(aipc_vals)-1); j0=int(v); j1=min(j0+1,len(aipc_vals)-1); dv=0 if j0==j1 else v-j0
            for kk in range(n):
                w=kk/(n-1)*(len(vhh_vals)-1); k0=int(w); k1=min(k0+1,len(vhh_vals)-1); dw=0 if k0==k1 else w-k0
                c000=cube[i0][j0][k0]; c100=cube[i1][j0][k0]; c010=cube[i0][j1][k0]; c110=cube[i1][j1][k0]
                c001=cube[i0][j0][k1]; c101=cube[i1][j0][k1]; c011=cube[i0][j1][k1]; c111=cube[i1][j1][k1]
                c00=c000*(1-du)+c100*du; c01=c001*(1-du)+c101*du
                c10=c010*(1-du)+c110*du; c11=c011*(1-du)+c111*du
                c0=c00*(1-dv)+c10*dv; c1=c01*(1-dv)+c11*dv
                c=c0*(1-dw)+c1*dw
                dens.append({
                    'x': ii/(n-1),
                    'y': jj/(n-1),
                    'z': kk/(n-1),
                    'Score': c
                })
    print(f'[info] dense {n}^3={n**3}')
    return dens

def project(x,y,z, rotY, rotX, perspective, zoom, W, H):
    cy=math.cos(rotY); sy=math.sin(rotY); cx=math.cos(rotX); sx=math.sin(rotX)
    xr=x*cy+z*sy; zr=-x*sy+z*cy
    yr=y*cx - zr*sx; zr2=y*sx + zr*cx
    if perspective:
        d=1/(1+zr2*0.8)
        return xr*zoom*d+W/2, -yr*zoom*d+H/2, zr2
    return xr*zoom+W/2, -yr*zoom+H/2, zr2

def render_frame(rows, dense, rotY, rotX, params):
    W,H=params.width, params.height
    img=Image.new('RGB',(W,H),params.bg_color)
    dr=ImageDraw.Draw(img)
    font_small=params.font_small; font_axis=params.font_axis
    scores=[r['Score'] for r in rows]; smin=min(scores); smax=max(scores)
    zoom=params.zoom; perspective=params.perspective
    origin=(-0.5,-0.5,-0.5)
    axes_info=[
        ((1,0,0),'vHH_to_cytokine_degradation_rate', params.dmin, params.dmax, (255,107,85)),
        ((0,1,0),'AIP_C_secretion_saturation',       params.amin, params.amax, (107,255,85)),
        ((0,0,1),'vHH_secretion_saturation',         params.vmin, params.vmax, (85,168,255)),
    ]
    ticks=params.ticks; w_axis=max(3,int(3*params.text_scale))
    for direction,label,pmin,pmax,color in axes_info:
        end=(origin[0]+direction[0],origin[1]+direction[1],origin[2]+direction[2])
        ox,oy,_=project(*origin,rotY,rotX,perspective,zoom,W,H)
        ex,ey,_=project(*end,rotY,rotX,perspective,zoom,W,H)
        dr.line((ox,oy,ex,ey),fill=color,width=w_axis)
        ah=int(10*params.text_scale); vx=ex-ox; vy=ey-oy; L=max(math.hypot(vx,vy),1e-6); ux,uy=vx/L,vy/L
        dr.polygon([(ex,ey),(ex-ux*ah - uy*ah*0.6, ey-uy*ah + ux*ah*0.6),(ex-ux*ah + uy*ah*0.6, ey-uy*ah - ux*ah*0.6)],fill=color)
        dr.text((ex+12*params.text_scale,ey),label,fill=color,font=font_axis,anchor='lm')
        for i in range(ticks):
            frac=i/(ticks-1 if ticks>1 else 1)
            tx=origin[0]+direction[0]*frac; ty=origin[1]+direction[1]*frac; tz=origin[2]+direction[2]*frac
            px_,py_,_ = project(tx,ty,tz,rotY,rotX,perspective,zoom,W,H)
            r_tick=int(3*params.text_scale)
            dr.ellipse((px_-r_tick,py_-r_tick,px_+r_tick,py_+r_tick),fill=(230,230,230))
            val=pmin+(pmax-pmin)*frac
            dr.text((px_,py_-16*params.text_scale),f"{val:.2f}",fill=(200,200,200),font=font_small,anchor='mm')
    if dense:
        for d in dense:
            x=d['x']-0.5; y=d['y']-0.5; z=d['z']-0.5
            px,py,_=project(x,y,z,rotY,rotX,perspective,zoom,W,H)
            t=(d['Score']-smin)/(smax-smin or 1); col=viridis(t)
            dr.ellipse((px-1,py-1,px+1,py+1),fill=col)
    # depth-sort original points
    deg_vals=params.deg_vals; aipc_vals=params.aipc_vals; vhh_vals=params.vhh_vals
    dmin,dmax=deg_vals[0],deg_vals[-1]; amin,amax=aipc_vals[0],aipc_vals[-1]; vmin,vmax=vhh_vals[0],vhh_vals[-1]
    pts=[]
    for r in rows:
        xn=(r['deg']-dmin)/(dmax-dmin or 1)-0.5; yn=(r['aipc']-amin)/(amax-amin or 1)-0.5; zn=(r['vhh']-vmin)/(vmax-vmin or 1)-0.5
        px,py,depth=project(xn,yn,zn,rotY,rotX,perspective,zoom,W,H)
        t=(r['Score']-smin)/(smax-smin or 1)
        pts.append((depth,px,py,t,r))
    pts.sort(key=lambda x:x[0])
    top_cut=int(len(rows)*params.highlight_top)
    top_ids=set(r['run_id'] for r in sorted(rows,key=lambda r:r['rank'])[:top_cut]) if top_cut>0 else set()
    base=5*params.text_scale; hi=7*params.text_scale
    for depth,px,py,t,r in pts:
        col=viridis(t)
        radius=(base if r['run_id'] not in top_ids else hi)*params.point_scale
        dr.ellipse((px-radius,py-radius,px+radius,py+radius),fill=col)
        if r['run_id'] in top_ids:
            off=params.text_scale*0.4
            dr.ellipse((px-radius-off,py-radius-off,px+radius+off,py+radius+off),outline=(255,255,255),width=max(1,int(params.text_scale)))
    if not params.no_legend:
        bar_w=int(160*params.text_scale); bar_h=int(14*params.text_scale); x0=20; y0=H-int(60*params.text_scale)
        for i in range(bar_w):
            tt=i/(bar_w-1); dr.line((x0+i,y0,x0+i,y0+bar_h),fill=viridis(tt))
        dr.rectangle((x0,y0,x0+bar_w,y0+bar_h),outline=(180,180,180))
        dr.text((x0,y0-12*params.text_scale),'Color = Score (low→high)',(220,220,220),font=font_small)
        dr.text((x0,y0+bar_h+2*params.text_scale),f"{smin:.3f}",(200,200,200),font=font_small)
        dr.text((x0+bar_w-60*params.text_scale,y0+bar_h+2*params.text_scale),f"{smax:.3f}",(200,200,200),font=font_small)
        dr.text((x0,20*params.text_scale),f"Axis={params.axis} Zoom={zoom} AngleScale={params.angle_scale}",(200,200,200),font=font_small)
        dr.text((x0,36*params.text_scale),'Top runs = white ring',(200,200,200),font=font_small)
    return img

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--csv',type=Path,default=Path('scored_runs.csv'))
    ap.add_argument('--output',type=Path,default=Path('rotation.gif'))
    ap.add_argument('--frames',type=int,default=120)
    ap.add_argument('--fps',type=int,default=24)
    ap.add_argument('--axis',type=str,default='y',choices=['x','y','z','oblique'])
    ap.add_argument('--sweep-angle',type=float,default=360.0,help='Base total yaw rotation (degrees).')
    ap.add_argument('--angle-scale',type=float,default=1.0,help='Multiply total rotation angle (0.1=1/10, 2=720°).')
    ap.add_argument('--start-elev',type=float,default=32.0,help='Starting elevation angle (degrees).')
    ap.add_argument('--zoom',type=float,default=480)
    ap.add_argument('--perspective',action='store_true')
    ap.add_argument('--dense',type=int,default=0)
    ap.add_argument('--highlight-top',type=float,default=0.05)
    ap.add_argument('--size',type=int,nargs=2,default=[900,900])
    ap.add_argument('--bg-color',type=str,default='#000000')
    ap.add_argument('--no-legend',action='store_true')
    ap.add_argument('--ticks',type=int,default=6,help='Ticks per axis (>=2).')
    ap.add_argument('--point-scale',type=float,default=1.0,help='Additional scale factor for point radii.')
    ap.add_argument('--font',type=Path,default=None,help='Optional TrueType font path.')
    ap.add_argument('--export-mp4',action='store_true',help='Also export MP4 via ffmpeg.')
    args=ap.parse_args()

    rows=read_rows(args.csv)
    if not rows:
        print('[ERROR] No rows read from',args.csv); sys.exit(1)
    dense=trilinear_dense(rows,args.dense) if args.dense>1 else []
    deg_vals=sorted({r['deg'] for r in rows}); aipc_vals=sorted({r['aipc'] for r in rows}); vhh_vals=sorted({r['vhh'] for r in rows})
    dmin,dmax=deg_vals[0],deg_vals[-1]; amin,amax=aipc_vals[0],aipc_vals[-1]; vmin,vmax=vhh_vals[0],vhh_vals[-1]

    class P: pass
    params=P()
    params.width,argsW=args.size[0],args.size[0]
    params.height=args.size[1]
    params.zoom=args.zoom; params.perspective=args.perspective; params.axis=args.axis
    params.highlight_top=args.highlight_top; params.no_legend=args.no_legend; params.bg_color=args.bg_color
    params.angle_scale=args.angle_scale; params.ticks=max(2,args.ticks)
    params.dmin,params.dmax=dmin,dmax; params.amin,params.amax=amin,amax; params.vmin,params.vmax=vmin,vmax
    params.deg_vals=deg_vals; params.aipc_vals=aipc_vals; params.vhh_vals=vhh_vals
    params.text_scale=3.0  # fixed 3x enlargement as requested
    params.point_scale=args.point_scale
    small_sz=int(12*params.text_scale); axis_sz=int(16*params.text_scale)
    try:
        if args.font and args.font.exists():
            params.font_small=ImageFont.truetype(str(args.font),small_sz)
            params.font_axis=ImageFont.truetype(str(args.font),axis_sz)
        else:
            params.font_small=ImageFont.load_default(); params.font_axis=ImageFont.load_default()
    except Exception:
        params.font_small=ImageFont.load_default(); params.font_axis=ImageFont.load_default()

    frames=[]; yaw0=0.0; elev0=math.radians(args.start_elev)
    total=math.radians(args.sweep_angle * args.angle_scale)
    print(f"[info] Generating {args.frames} frames, total yaw {math.degrees(total):.1f}° axis={args.axis}")
    for f in range(args.frames):
        t=f/(args.frames-1 if args.frames>1 else 1)
        if args.axis=='y':
            rotY=yaw0+total*t; rotX=elev0
        elif args.axis=='x':
            rotY=yaw0; rotX=elev0+total*t
        elif args.axis=='z':
            rotY=yaw0+total*t; rotX=elev0  # treat like yaw for simplicity
        else: # oblique
            rotY=yaw0+total*t; rotX=elev0+total*0.6*t
        img=render_frame(rows,dense,rotY,rotX,params)
        frames.append(img)
        if (f+1)%10==0: print(f'  frame {f+1}/{args.frames}')

    duration=int(1000/args.fps)
    frames[0].save(args.output,save_all=True,append_images=frames[1:],loop=0,duration=duration,disposal=2)
    print('[done] GIF written:',args.output)

    if args.export_mp4:
        if shutil.which('ffmpeg'):
            mp4=args.output.with_suffix('.mp4')
            cmd=['ffmpeg','-y','-i',str(args.output),'-movflags','faststart','-pix_fmt','yuv420p',str(mp4)]
            print('[info] Exporting MP4 ...')
            try:
                subprocess.run(cmd,check=True)
                print('[done] MP4 written:',mp4)
            except Exception as e:
                print('[warn] MP4 export failed:',e)
        else:
            print('[warn] ffmpeg not found; skipping MP4 export')
    else:
        print('Hint: add --export-mp4 to also create an MP4 file')

    print('Done.')

if __name__=='__main__':
    main()
