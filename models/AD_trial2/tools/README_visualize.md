Visualize PhysiCell .mat outputs

This small tool helps you turn PhysiCell/BioFVM .mat microenvironment outputs into:

- Per-timepoint heatmap PNGs (single z-slice or max-projection)
- MP4 videos (one per substrate)
- Voxel-averaged concentration vs time plots (multi-substrate)

Requirements

- Python 3.8+
- scipy
- numpy
- matplotlib

Install (PowerShell example)

```powershell
python -m pip install --upgrade pip; python -m pip install scipy numpy matplotlib
```

Basic usage

```powershell
# generate PNG heatmaps (max-projection) and a time-series plot for substrates 0 and 1
python tools\visualize_physicell_mat.py .\output --substrates 0 1 --grid-shape 50 50 10 --outdir viz_out --plot

# also produce MP4 videos for those substrates
python tools\visualize_physicell_mat.py .\output --substrates 0 1 --grid-shape 50 50 10 --outdir viz_out --plot --video --fps 8
```

Notes

- You must provide `--grid-shape NX NY NZ` when the .mat files contain voxel vectors (typical). The script will reshape the 1D concentration vector into an (NX,NY,NZ) grid using Fortran ('F') ordering which matches many PhysiCell outputs. If your grid uses C-ordering, swap to order='C' in `make_heatmap_grid`.
- Substrates can be provided by index (0-based). If your .mat files include substrate names in a field called `substrate_names`, the script will try to read them (but this is not guaranteed for all output formats).
- If the script cannot find the substrate arrays automatically, open one of your `.mat` files in Python and inspect keys to find where the microenvironment is stored.

Troubleshooting

- "No .mat files found": point `mat_folder` at the folder that contains the `.mat` files (not the parent folder).
- "Could not reshape to grid": verify your `--grid-shape` matches the simulation mesh.
- If you want a different colormap or vmin/vmax scaling, edit `save_heatmap` in the script.

Contact

If your .mat files use a different structure, send one (or a redacted example) and I can extend the extractor to match it.

Additional examples

- Create an MP4 from PNGs using `ffmpeg` (PowerShell):

```powershell
ffmpeg -framerate 8 -pattern_type glob -i "viz_out/sub_0_t*.png" -c:v libx264 -pix_fmt yuv420p viz_out/sub_0.mp4
```

- Example: run the tool for substrates 2 and 3 with a different grid, saving results to `viz_results`:

```powershell
python tools\visualize_physicell_mat.py .\output --substrates 2 3 --grid-shape 100 100 5 --outdir viz_results --plot --video
```
