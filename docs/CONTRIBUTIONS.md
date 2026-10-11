# Personal contribution and provenance

**Yuling Zhou was the team's sole modeling member**, as stated by the portfolio owner. This repository presents the modeling work carried out for SJTU-BioX-Shanghai during the 2025 iGEM project.

| Contribution | Evidence to inspect |
| --- | --- |
| Molecular ODE formulation, circuit dynamics and sensitivity interpretation | [Molecular wiki chapter](wiki/model-part1.md): reaction equations, parameter sources, trajectories and sensitivity plots |
| Microbial interaction and host-response rules | `models/*/config/cell_rules.csv`; [spatial modeling description](RESEARCH.md) |
| Custom substrate dynamics and biofilm/barrier capture | `models/*/custom_modules/custom.cpp`, especially `update_substrate_reactions` and EPS/ECM handling |
| 2D occurrence, recurrence and treatment scenarios | `AD_trial1`, `AD_trial2`, `Recurrence`, `Treatment`: initial cell CSVs and active settings |
| 3D tissue layout and colony placement | `AD_3D/tools/` and `Recurrence_3D/config/generate_engr_epidermis*.py`; final supplied `cells.csv` files |
| Population and substrate analysis, bivariate Moran's I | `AD_trial2/tools/extract_population_timeseries.py` and `visualize_physicell_mat.py` |
| Parameter exploration and design trade-offs | `Recurrence/parameter_sweep.py`, `score_runs.py` and 1,000 archived run folders |

## Contribution boundaries

The role statement describes responsibility for the modeling work. It does not claim authorship of every line of code in the repository. PhysiCell/BioFVM, Hill-response rule infrastructure, and other upstream software retain their own authors and license notices. Some original helper scripts contain assistant-generated attribution; those headers are retained. New portfolio assembly, portable execution and verification scripts were prepared with coding-assistant support in October 2026.

Wet-lab assays, the broader therapeutic concept and team project outcomes are collaborative team contributions. Experimental measurements cited in the drafts are distinguished from simulations. No award, publication or clinical efficacy claim has been added to this portfolio.

## Original sources

- [Public iGEM GitLab project](https://gitlab.igem.org/linnazhou/model-tools-for-2025-sjtu-biox-shanghai), source commit `6d03dda30d6ff2c6e825f2a66b8d3ef4b4a549f5`, dated 2025-10-08.
- Local `wiki_description/`: six original Word documents, copied byte-for-byte to [originals](originals).
- Local PhysiCell 1.14.2: core, modules, BioFVM and license files, retained separately under `vendor/PhysiCell`.

[source-manifest.json](source-manifest.json) records source locations and file hashes, plus excluded backups and runtime artifacts. The original GitLab repository remains the reference for historical commits. This curated repository is a documented portfolio export rather than a full GitLab history mirror.
