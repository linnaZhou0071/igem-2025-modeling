#!/usr/bin/env python3
"""
visualize_physicell_mat.py

Tool to visualize PhysiCell/BioFVM .mat microenvironment outputs.

Features:
- Read a folder of .mat files produced by PhysiCell (one per output timepoint)
- Produce heatmap PNGs for selected substrates at each timepoint (single z-slice or max-projection)
- Produce an MP4 video per substrate
- Produce voxel-averaged concentration vs time plots for selected substrates (multiple on one plot)

Usage examples are in the accompanying README.
"""
import argparse
import os
import sys
import glob
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import animation
from scipy import sparse
import random
import csv  # for optional --series-csv output


def find_mat_files(folder):
    # Find candidate .mat files and inspect them to ensure they contain
    # microenvironment-style arrays (compact layout or M/microenvironment vars).
    all_files = sorted(glob.glob(os.path.join(folder, "*.mat")))
    candidates = []
    for fp in all_files:
        base = os.path.basename(fp).lower()
        # quick skip for obvious cell-only files
        if 'cells.mat' in base or base.startswith('output') and base.endswith('cells.mat'):
            continue
        # try to inspect variable info without fully loading large arrays
        try:
            from scipy.io import whosmat
            info = whosmat(fp)
            # info is list of (name, shape, dtype)
            is_me = False
            for name, shape, dtype in info:
                lname = name.lower()
                # variable name hints
                if 'micro' in lname or name in ('M', 'M0', 'microenvironment', 'multiscale_microenvironment') or 'substrate' in lname or 'density' in lname:
                    is_me = True
                    break
                # shape hints: compact layout often has rows >= 4 (3 coords + substrates)
                if shape is not None and len(shape) == 2:
                    r, c = shape
                    if (r >= 4) or (c >= 4):
                        is_me = True
                        break
            if is_me:
                candidates.append(fp)
        except Exception:
            # if whosmat fails, conservatively include the file by name pattern
            if 'microenvironment' in base or 'microenv' in base:
                candidates.append(fp)

    if not candidates:
        print('No microenvironment-like .mat files detected in', folder)
    return sorted(candidates)


def load_mat_file(path):
    data = scipy.io.loadmat(path)
    return data


def list_substrates(mat):
    # Common PhysiCell/Microenvironment variable naming
    # Try multiple keys to be robust
    for key in ["microenv", "microenvironment", "M", "microenvironment_variables"]:
        if key in mat:
            m = mat[key]
            break
    # fallback: inspect keys
    keys = list(mat.keys())
    # Heuristic: look for 'substrate_names' or 'names'
    if "substrate_names" in mat:
        names = [str(x[0]) for x in mat["substrate_names"].ravel()]
        return names
    for k in keys:
        if "name" in k.lower() and isinstance(mat[k], np.ndarray):
            try:
                vals = mat[k]
                names = [str(x[0]) for x in vals.ravel()]
                return names
            except Exception:
                pass
    # Last resort: user must provide indices
    return None


def extract_substrate_grid(mat, substrate_index, grid_shape=None):
    """
    Robustly search a loaded .mat dict for a substrate array and return a 1D
    array for the requested substrate_index. When grid_shape is provided we
    use the voxel count to identify arrays of the expected length.
    """
    nvox = None
    if grid_shape is not None:
        try:
            nvox = int(grid_shape[0]) * int(grid_shape[1]) * int(grid_shape[2])
        except Exception:
            nvox = None

    def inspect_obj(obj):
        # Handle MATLAB structs loaded as object arrays
        if isinstance(obj, np.ndarray) and obj.dtype == object:
            try:
                # many structs are (1,1) object arrays
                inner = obj[0, 0]
                # if inner is a numpy structured array with fields
                if hasattr(inner, 'dtype') and getattr(inner.dtype, 'names', None):
                    for fname in inner.dtype.names:
                        try:
                            f = inner[fname]
                        except Exception:
                            continue
                        res = inspect_obj(f)
                        if res is not None:
                            return res
                # otherwise, recurse into the inner object
                return inspect_obj(inner)
            except Exception:
                pass

        try:
            a = np.array(obj)
        except Exception:
            return None

        if not isinstance(a, np.ndarray):
            return None

        # Handle the compact MATLAB microenvironment layout described in docs:
        # one column per voxel, rows = 3 + n_substrates (rows 0..2 = x,y,z)
        if a.ndim == 2:
            r, c = a.shape
            # case: rows = 3 + n_substrates (shape = (nrows, nvoxels))
            if r >= 4:
                nsubstrates_here = r - 3
                if substrate_index < nsubstrates_here:
                    # return the row corresponding to this substrate
                    return a[3 + substrate_index, :]
            # case: transposed (shape = (nvoxels, nrows))
            if c >= 4:
                nsubstrates_here = c - 3
                if substrate_index < nsubstrates_here:
                    return a[:, 3 + substrate_index]

        # 1D arrays: possibly a single-substrate voxel vector
        if a.ndim == 1:
            if nvox is None or a.size == nvox:
                if substrate_index == 0:
                    return a

        return None

    # search top-level entries
    for key in mat:
        if key.startswith('__'):
            continue
        res = inspect_obj(mat[key])
        if res is not None:
            return np.array(res).ravel()

    raise ValueError('Could not locate substrate array in .mat file; please provide substrate indices or a file with names')


def voxel_average_from_array(arr):
    return np.mean(arr)


def plot_time_series(times, series_dict, outpath, title=None):
    plt.figure(figsize=(8,5))
    for name, vals in series_dict.items():
        plt.plot(times, vals, label=name)
    plt.xlabel('time')
    plt.ylabel('voxel-averaged concentration')
    plt.legend()
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


def make_heatmap_grid(data_1d, grid_shape, slice_z=None):
    # data_1d is length n_voxels; grid_shape is (nx, ny, nz)
    nx, ny, nz = grid_shape
    expected = int(nx) * int(ny) * int(nz)
    actual = int(np.array(data_1d).ravel().size)
    if actual != expected:
        # try 2D fallback: if actual == nx*ny and nz is 1 or the simulation is 2D
        if actual == int(nx) * int(ny):
            arr3 = np.zeros((nx, ny, max(1, nz)))
            arr3[:, :, 0] = np.array(data_1d).ravel().reshape((nx, ny), order='F')
        else:
            raise ValueError(f'voxel count mismatch: expected {expected}, got {actual}')
    else:
        arr3 = data_1d.reshape((nx, ny, nz), order='F')
    if slice_z is None:
        # max projection
        img = np.max(arr3, axis=2)
    else:
        img = arr3[:, :, slice_z]
    return img


def save_heatmap(img, outpath, cmap='viridis', vmin=None, vmax=None):
    plt.figure(figsize=(6,5))
    plt.imshow(img.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Visualize PhysiCell .mat microenvironment outputs')
    parser.add_argument('mat_folder', help='Folder containing .mat files')
    parser.add_argument('--substrates', '-s', nargs='+', required=True,
                        help='Substrate indices (0-based) or names (if available)')
    parser.add_argument('--grid-shape', '-g', nargs=3, type=int, metavar=('NX','NY','NZ'),
                        help='Grid shape (required to reconstruct 3D from voxel vectors)')
    parser.add_argument('--slice-z', type=int, default=None, help='Z slice index to plot (default: maxprojection)')
    parser.add_argument('--outdir', '-o', default='viz_out', help='Output folder')
    parser.add_argument('--video', action='store_true', help='Produce MP4 video per substrate')
    parser.add_argument('--plot', action='store_true', help='Produce averaged time-series plot')
    parser.add_argument('--fps', type=int, default=6, help='Video frames per second')
    parser.add_argument('--one-based', action='store_true', help='Treat numeric substrate indices as 1-based (convert to 0-based)')
    parser.add_argument('--no-images', action='store_true', help='Do not save per-timestep heatmap images; only compute averages and plots')
    parser.add_argument('--timeplot', nargs='*', help='List of substrates (indices or names) to include on the time-series plot (defaults to all substrates)')
    # --- Bivariate Moran's I extensions ---
    parser.add_argument('--bivar-moran-pairs', '-bm', nargs='*',
                        help='Pairs of substrates for global bivariate Moran\'s I over time. Format examples: 4,10  1:2  3-7')
    parser.add_argument('--bivar-moran-norm', default='xy', choices=['xy','x'],
                        help="Normalization: 'xy' uses sqrt(Var(x)Var(y)) giving a correlation-like coefficient; 'x' uses Var(x) only (asymmetric, like classic univariate Moran form).")
    parser.add_argument('--bivar-moran-perm', type=int, default=0,
                        help='If >0, run this many permutations (randomly shuffling first substrate values) per timepoint to estimate p-values.')
    parser.add_argument('--bivar-moran-output', default='bivar_moran',
                        help='Prefix (within outdir) for bivariate Moran outputs (CSV + PNG).')
    # --- Validation / diagnostics (added) ---
    parser.add_argument('--validate', action='store_true',
                        help='Run extra sanity checks (orientation, voxel count, NaNs, monotonic time) and print summary.')
    parser.add_argument('--stats-csv', default=None,
                        help='If set, write per-timepoint summary stats (mean/std/min/max/percentiles) for each substrate.')
    parser.add_argument('--track-voxels', nargs='+', type=int,
                        help='Linear voxel indices (0-based) to track raw concentrations over time (requires --grid-shape).')
    parser.add_argument('--track-output', default='tracked_voxels.csv',
                        help='Output CSV filename for --track-voxels data (default: tracked_voxels.csv).')
    parser.add_argument('--series-csv', default=None,
                        help='Write a wide CSV of voxel-averaged concentrations over time (columns: time, one column per requested substrate). No images needed; combine with --no-images and omit --plot/--video for pure numeric export.')
    parser.add_argument('--export-voxels', action='store_true',
                        help='Export per-timepoint, per-substrate voxel concentration CSVs. 2D (nz=1): matrix form (nx rows x ny cols). 3D (nz>1): long form with columns x,y,z,value.')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    mat_files = find_mat_files(args.mat_folder)
    if not mat_files:
        print('No .mat files found in', args.mat_folder)
        sys.exit(1)

    # Try to load substrate names once (even if we'll skip time-zero later)
    substrate_names = None
    try:
        _first_mat = load_mat_file(mat_files[0])
        substrate_names = list_substrates(_first_mat)
    except Exception:
        substrate_names = None

    # parse substrate selection: convert names to indices when possible
    substrates = []
    try:
        substrates = [int(x) for x in args.substrates]
        if args.one_based:
            substrates = [int(x) - 1 for x in substrates]
    except Exception:
        substrates = args.substrates

    times = []
    series = {str(s): [] for s in substrates}
    # For bivariate Moran: store full field per substrate per timepoint (flattened voxel vectors)
    store_fields = (args.bivar_moran_pairs is not None and len(args.bivar_moran_pairs) > 0)
    fields_time_series = {str(s): [] for s in substrates} if store_fields else None

    # For video: store images per substrate
    frames_per_substrate = {str(s): [] for s in substrates}

    # Validation / tracking state containers
    nvox_expected = None
    if args.grid_shape is not None:
        try:
            nvox_expected = int(args.grid_shape[0]) * int(args.grid_shape[1]) * int(args.grid_shape[2])
        except Exception:
            nvox_expected = None
    stats_rows = []  # (time, substrate, nvox, nNaN, mean, std, min, p25, median, p75, max)
    track_indices = args.track_voxels if args.track_voxels else []
    if track_indices and nvox_expected is None:
        print('[WARN] --track-voxels provided but --grid-shape missing; cannot fully validate indices.')
    tracked_series = {str(s): {idx: [] for idx in track_indices} for s in substrates}
    validation = {
        'errors': 0,
        'warnings': 0,
        'first_timepoint_checked': False,
        'nan_substrates': set(),
    }

    for path in mat_files:
        mat = load_mat_file(path)
        # try to parse time from filename (common: output_t0000.mat)
        basename = os.path.basename(path)
        import re
        m = re.search(r't(\d+)', basename)
        if m:
            t = int(m.group(1))
        else:
            t = len(times)
        # skip initial / time-zero files so they do not generate PNGs or affect averages/videos
        if 'initial' in basename.lower() or t == 0:
            print('Skipping initial/time-zero file', basename)
            continue
        times.append(t)

        # temporary dict to store this timestep's flattened arrays (for speed / avoid recomputation)
        flat_cache = {}
        for s in substrates:
            key = str(s)
            try:
                arr = extract_substrate_grid(mat, int(s), grid_shape=args.grid_shape)
            except Exception as e:
                # print available keys and brief shapes to help debugging
                try:
                    keys = [k for k in mat.keys() if not k.startswith('__')]
                    info = []
                    for k in keys:
                        try:
                            a = np.array(mat[k])
                            info.append(f"{k}:{getattr(a,'shape',None)}")
                        except Exception:
                            info.append(f"{k}:<non-array>")
                    print('Failed to extract substrate', s, 'from', path, ':', e, '; keys:', info)
                except Exception:
                    print('Failed to extract substrate', s, 'from', path, ':', e)
                arr = None
            if arr is None:
                series[key].append(np.nan)
                frames_per_substrate[key].append(None)
                if args.stats_csv is not None:
                    stats_rows.append([t, key, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
                continue

            # arr may be 1D vector of length n_voxels
            if args.grid_shape:
                try:
                    img = make_heatmap_grid(np.array(arr).ravel(), tuple(args.grid_shape), slice_z=args.slice_z)
                except Exception as e:
                    print('Could not reshape to grid for', path, ':', e)
                    img = None
            else:
                img = None

            # (NEW) voxel CSV export
            if args.export_voxels and args.grid_shape:
                nx, ny, nz = args.grid_shape
                try:
                    arr3 = np.array(arr).ravel().reshape((nx, ny, nz), order='F')
                    # Decide filename pattern; include substrate name if available
                    disp_name = None
                    if substrate_names is not None:
                        try:
                            ik = int(key)
                            if 0 <= ik < len(substrate_names):
                                disp_name = substrate_names[ik]
                        except Exception:
                            pass
                    if disp_name is not None:
                        # sanitize name
                        safe_name = ''.join(ch if ch.isalnum() or ch in ('_','-') else '_' for ch in disp_name)
                        voxel_fname = os.path.join(args.outdir, f'vox_t{t:04d}_sub{key}_{safe_name}.csv')
                    else:
                        voxel_fname = os.path.join(args.outdir, f'vox_t{t:04d}_sub{key}.csv')
                    if nz == 1:
                        # matrix form
                        with open(voxel_fname, 'w', newline='') as vf:
                            # optional header
                            vf.write(f'# substrate={key} time={t} shape=({nx},{ny},{nz}) format=matrix\n')
                            writer = csv.writer(vf)
                            for x_ in range(nx):
                                row_vals = []
                                for y_ in range(ny):
                                    v = arr3[x_, y_, 0]
                                    try:
                                        v = float(np.asarray(v).ravel()[0])
                                    except Exception:
                                        v = np.nan
                                    row_vals.append(v)
                                writer.writerow([f'{v:.6g}' if np.isfinite(v) else '' for v in row_vals])
                    else:
                        # long form
                        with open(voxel_fname, 'w', newline='') as vf:
                            vf.write(f'# substrate={key} time={t} shape=({nx},{ny},{nz}) format=long(x,y,z,value)\n')
                            writer = csv.writer(vf)
                            writer.writerow(['x','y','z','value'])
                            for x_ in range(nx):
                                for y_ in range(ny):
                                    for z_ in range(nz):
                                        v = arr3[x_, y_, z_]
                                        try:
                                            v = float(np.asarray(v).ravel()[0])
                                        except Exception:
                                            v = np.nan
                                        writer.writerow([x_, y_, z_, f'{v:.6g}' if np.isfinite(v) else ''])
                    # simple progress output for large batches
                    print(f'[voxels] wrote {voxel_fname}')
                except Exception as e:
                    print(f'[WARN] failed voxel export for substrate {key} time {t}:', e)

            if img is not None:
                if not args.no_images:
                    out_png = os.path.join(args.outdir, f'sub_{s}_t{t:04d}.png')
                    save_heatmap(img, out_png)
                    frames_per_substrate[key].append(out_png)
                else:
                    frames_per_substrate[key].append(None)
            else:
                frames_per_substrate[key].append(None)

            # voxel-average
            try:
                avg = voxel_average_from_array(np.array(arr).ravel())
            except Exception:
                avg = np.nan
            series[key].append(float(avg))
            if store_fields:
                if arr is not None:
                    flat_cache[key] = np.array(arr).ravel().astype(float)
                else:
                    flat_cache[key] = None

            # --- Validation & stats (added) ---
            if args.validate:
                # Check voxel count consistency
                if nvox_expected is not None and arr.size != nvox_expected:
                    print(f'[ERROR] time {t} substrate {key}: vector length {arr.size} != expected {nvox_expected}.')
                    validation['errors'] += 1
                # Very small length heuristic warning (likely orientation issue if >4 substrates expected)
                if nvox_expected is None and arr.size < 10:
                    print(f'[WARN] time {t} substrate {key}: suspiciously small vector length {arr.size}.')
                    validation['warnings'] += 1
                # NaN check
                nan_ct = int(np.isnan(arr).sum())
                if nan_ct > 0:
                    print(f'[WARN] time {t} substrate {key}: {nan_ct} NaN voxels.')
                    validation['warnings'] += 1
                    validation['nan_substrates'].add(key)
            # Summary stats
            if args.stats_csv is not None:
                a = np.array(arr).astype(float)
                nan_mask = np.isnan(a)
                n_nan = int(nan_mask.sum())
                if n_nan == a.size:
                    stats_rows.append([t, key, a.size, n_nan] + [np.nan]*7)
                else:
                    # use nan* functions for robustness
                    valid = a[~nan_mask]
                    p25, med, p75 = np.percentile(valid, [25, 50, 75])
                    stats_rows.append([t, key, a.size, n_nan,
                                       float(np.nanmean(a)), float(np.nanstd(a)),
                                       float(np.nanmin(a)), float(p25), float(med), float(p75), float(np.nanmax(a))])
            # Tracking specific voxels
            if track_indices:
                for idx in track_indices:
                    if idx < 0 or idx >= arr.size:
                        tracked_series[key][idx].append(np.nan)
                        if args.validate:
                            print(f'[WARN] time {t} substrate {key}: track index {idx} out of bounds (size {arr.size}).')
                            validation['warnings'] += 1
                    else:
                        tracked_series[key][idx].append(float(arr[idx]))
        if store_fields:
            # ensure each substrate list stays aligned with 'times'
            for k in fields_time_series.keys():
                fields_time_series[k].append(flat_cache.get(k, None))

    # write time-series plot
    if args.plot:
        # determine which series to plot
        if args.timeplot:
            # parse requested plot keys (may be indices or names)
            plot_keys = []
            for x in args.timeplot:
                try:
                    ix = int(x)
                    if args.one_based:
                        ix = ix - 1
                    plot_keys.append(str(ix))
                except Exception:
                    plot_keys.append(str(x))
        else:
            plot_keys = [str(s) for s in substrates]

        series_to_plot = {}
        for k in plot_keys:
            if k in series:
                series_to_plot[k] = series[k]
            else:
                print('Warning: requested plot substrate', k, "not found in computed series")

        out_ts = os.path.join(args.outdir, 'time_series.png')
        plot_time_series(times, series_to_plot, out_ts)
        print('Wrote', out_ts)

    # NEW: plain wide time-series CSV of voxel-averaged concentrations
    if args.series_csv:
        out_csv_path = os.path.join(args.outdir, args.series_csv)
        try:
            with open(out_csv_path, 'w', newline='') as f:
                w = csv.writer(f)
                header = ['time'] + list(series.keys())
                w.writerow(header)
                # Each series list aligns with times[] (skip any initial file we skipped above)
                for i, t in enumerate(times):
                    row = [t]
                    for k in series.keys():
                        vals = series[k]
                        row.append(vals[i] if i < len(vals) else '')
                    w.writerow(row)
            print('[INFO] Wrote series CSV ->', out_csv_path)
        except Exception as e:
            print('[ERROR] failed writing series CSV:', e)

    # write videos
    if args.video:
        for s in substrates:
            key = str(s)
            frames = frames_per_substrate[key]
            images = [f for f in frames if f is not None]
            if not images:
                print('No images for substrate', s, '; skipping video')
                continue
            # create animation
            fig = plt.figure()
            im = plt.imshow(plt.imread(images[0]).astype(float)/255.0)
            plt.axis('off')

            def update_frame(i):
                img = plt.imread(images[i])
                im.set_data(img.astype(float)/255.0)
                return (im,)

            ani = animation.FuncAnimation(fig, update_frame, frames=len(images), interval=1000/args.fps)
            out_mp4 = os.path.join(args.outdir, f'sub_{s}.mp4')
            ani.save(out_mp4, fps=args.fps, dpi=150)
            plt.close(fig)
            print('Wrote', out_mp4)

    # ---------------------------------------------------------------
    # Bivariate Moran's I (global) over time for specified pairs
    # ---------------------------------------------------------------
    if store_fields:
        # parse pairs into list of (a,b) substrate indices (as strings consistent with keys)
        def parse_pair(token):
            for sep in [',', ':', '-', '/']:
                if sep in token:
                    a,b = token.split(sep,1)
                    return a.strip(), b.strip()
            raise ValueError(f'Invalid pair format: {token}')

        pair_list = []
        for tok in args.bivar_moran_pairs:
            try:
                a,b = parse_pair(tok)
                # convert to indices (respect one-based if needed)
                try:
                    ai = int(a); bi = int(b)
                    if args.one_based:
                        ai -= 1; bi -= 1
                    pair_list.append( (str(ai), str(bi), tok) )
                except Exception:
                    pair_list.append( (a,b,tok) )
            except Exception as e:
                print('Skipping pair token', tok, ':', e)
        if not pair_list:
            print('No valid bivariate Moran pairs parsed.')
        else:
            if not args.grid_shape:
                print('Grid shape is required for bivariate Moran computation. Aborting that part.')
            else:
                nx, ny, nz = args.grid_shape
                nvox = int(nx)*int(ny)*int(nz)
                # Build sparse row-standardized 4/6-neighborhood weight matrix W
                def build_W(nx, ny, nz):
                    rows = []
                    cols = []
                    data = []
                    def idx(x,y,z):
                        return x + nx*(y + ny*z)
                    use_3d = (nz > 1)
                    neigh_2d = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0)]
                    neigh_3d = neigh_2d + [(0,0,1),(0,0,-1)]
                    neigh_list = neigh_3d if use_3d else neigh_2d
                    for x in range(nx):
                        for y in range(ny):
                            for z in range(nz):
                                i = idx(x,y,z)
                                deg = 0
                                for dx,dy,dz in neigh_list:
                                    xx,yy,zz = x+dx, y+dy, z+dz
                                    if 0 <= xx < nx and 0 <= yy < ny and 0 <= zz < nz:
                                        j = idx(xx,yy,zz)
                                        rows.append(i); cols.append(j); data.append(1.0)
                                        deg += 1
                    # row-standardization
                    W = sparse.coo_matrix((data,(rows,cols)), shape=(nvox,nvox)).tocsr()
                    row_sums = np.array(W.sum(axis=1)).ravel()
                    row_sums[row_sums==0] = 1.0
                    inv = 1.0/row_sums
                    W = sparse.diags(inv) @ W
                    return W
                W = build_W(nx, ny, nz)
                S0 = nvox  # after row-standardization

                def bivar_moran(x, y, W, S0, norm_mode='xy'):
                    # x, y are 1D arrays length nvox, may contain NaN
                    mask = ~(np.isnan(x) | np.isnan(y))
                    x = x[mask]; y = y[mask]
                    if x.size == 0:
                        return np.nan
                    # restrict W
                    if mask.sum() != W.shape[0]:
                        # create submatrix (this can be expensive; assume usually no NaNs)
                        W_sub = W[mask][:,mask]
                    else:
                        W_sub = W
                    x_mean = x.mean(); y_mean = y.mean()
                    x_c = x - x_mean; y_c = y - y_mean
                    # numerator: x_c^T W y_c
                    num = x_c @ (W_sub @ y_c)
                    if norm_mode == 'xy':
                        denom = np.sqrt( (x_c @ x_c) * (y_c @ y_c) )
                    else:  # 'x'
                        denom = (x_c @ x_c)
                    if denom == 0:
                        return np.nan
                    return (x.size / S0) * (num / denom)

                # Prepare output containers
                moran_results = {p[2]: [] for p in pair_list}  # key is original token
                perm_pvalues = {p[2]: [] for p in pair_list} if args.bivar_moran_perm > 0 else None

                for ti, t in enumerate(times):
                    for (a,b,token) in pair_list:
                        xa = fields_time_series.get(a, None)
                        yb = fields_time_series.get(b, None)
                        if xa is None or yb is None:
                            moran_results[token].append(np.nan)
                            if perm_pvalues: perm_pvalues[token].append(np.nan)
                            continue
                        x_vec = xa[ti]; y_vec = yb[ti]
                        if x_vec is None or y_vec is None:
                            moran_results[token].append(np.nan)
                            if perm_pvalues: perm_pvalues[token].append(np.nan)
                            continue
                        val = bivar_moran(x_vec, y_vec, W, S0, norm_mode=args.bivar_moran_norm)
                        moran_results[token].append(val)
                        # permutations for p-value (shuffle x only)
                        if perm_pvalues:
                            if np.isnan(val):
                                perm_pvalues[token].append(np.nan)
                            else:
                                count_ge = 0
                                valid_mask = ~(np.isnan(x_vec) | np.isnan(y_vec))
                                xv = x_vec[valid_mask]; yv = y_vec[valid_mask]
                                # Precompute basics for speed
                                W_sub = W if valid_mask.sum()==W.shape[0] else W[valid_mask][:,valid_mask]
                                x_mean = xv.mean(); y_mean = yv.mean()
                                x_c = xv - x_mean; y_c = yv - y_mean
                                base_denom_xy = np.sqrt( (x_c @ x_c) * (y_c @ y_c) ) if args.bivar_moran_norm=='xy' else (x_c @ x_c)
                                for r in range(args.bivar_moran_perm):
                                    perm = np.random.permutation(xv)
                                    perm_c = perm - perm.mean()  # re-center after shuffle
                                    num_r = perm_c @ (W_sub @ y_c)
                                    if args.bivar_moran_norm=='xy':
                                        denom_r = np.sqrt( (perm_c @ perm_c) * (y_c @ y_c) )
                                    else:
                                        denom_r = (perm_c @ perm_c)
                                    if denom_r == 0:
                                        continue
                                    ir = (xv.size / S0) * (num_r / denom_r)
                                    if ir >= val:
                                        count_ge += 1
                                pval = (count_ge + 1) / (args.bivar_moran_perm + 1)
                                perm_pvalues[token].append(pval)

                # Save CSV
                # csv already imported at module level
                csv_path = os.path.join(args.outdir, args.bivar_moran_output + '.csv')
                with open(csv_path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    header = ['time'] + [f'I({tok})' for tok in moran_results.keys()]
                    if perm_pvalues:
                        header += [f'p({tok})' for tok in perm_pvalues.keys()]
                    writer.writerow(header)
                    for i, t in enumerate(times):
                        row = [t]
                        for tok in moran_results.keys():
                            row.append(moran_results[tok][i])
                        if perm_pvalues:
                            for tok in perm_pvalues.keys():
                                row.append(perm_pvalues[tok][i])
                        writer.writerow(row)
                print('Wrote', csv_path)

                # Plot
                plt.figure(figsize=(8,5))
                for tok, vals in moran_results.items():
                    plt.plot(times, vals, label=f'I({tok})')
                plt.xlabel('time')
                if args.bivar_moran_norm=='xy':
                    plt.ylabel('Bivariate Moran I (correlation-like)')
                else:
                    plt.ylabel('Bivariate Moran I (variance-x normalized)')
                plt.legend()
                plt.title('Bivariate Moran\'s I over time')
                plt.tight_layout()
                png_path = os.path.join(args.outdir, args.bivar_moran_output + '.png')
                plt.savefig(png_path, dpi=150)
                plt.close()
                print('Wrote', png_path)

    # ---------------------------
    # Write stats / tracking output (added)
    # ---------------------------
    if args.stats_csv is not None:
        # csv already imported at module level
        header = ['time','substrate','nvox','nNaN','mean','std','min','p25','median','p75','max']
        try:
            with open(os.path.join(args.outdir, args.stats_csv), 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(header)
                for row in stats_rows:
                    w.writerow(row)
            print('[INFO] Wrote stats CSV ->', os.path.join(args.outdir, args.stats_csv))
        except Exception as e:
            print('[ERROR] failed writing stats CSV:', e)

    if track_indices:
        # csv already imported at module level
        # Build wide table: columns = time, then for each substrate+voxel a column name sub{substrate}_v{idx}
        col_names = ['time']
        for s in substrates:
            sk = str(s)
            for idx in track_indices:
                col_names.append(f'sub{sk}_v{idx}')
        try:
            with open(os.path.join(args.outdir, args.track_output), 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(col_names)
                for ti, t in enumerate(times):
                    row = [t]
                    for s in substrates:
                        sk = str(s)
                        for idx in track_indices:
                            series_list = tracked_series[sk][idx]
                            if ti < len(series_list):
                                row.append(series_list[ti])
                            else:
                                row.append(np.nan)
                    w.writerow(row)
            print('[INFO] Wrote tracked voxels CSV ->', os.path.join(args.outdir, args.track_output))
        except Exception as e:
            print('[ERROR] failed writing tracked voxels CSV:', e)

    # ---------------------------
    # Validation summary (added)
    # ---------------------------
    if args.validate:
        # Check time monotonicity
        if any(t2 <= t1 for t1, t2 in zip(times, times[1:])):
            print('[ERROR] Non-monotonic or duplicate time values detected.')
            validation['errors'] += 1
        print('==== Validation Summary ====')
        print('Time points processed :', len(times))
        if nvox_expected is not None:
            print('Expected voxels       :', nvox_expected)
        print('Warnings              :', validation['warnings'])
        print('Errors                :', validation['errors'])
        if validation['nan_substrates']:
            print('Substrates with NaNs  :', ', '.join(sorted(validation['nan_substrates'])))
        if validation['errors'] == 0:
            print('Overall               : PASS')
        else:
            print('Overall               : FAIL')


if __name__ == '__main__':
    main()
