# Editorial and evidence notes

The main README and research overview are curated portfolio descriptions. The locally retained Word files preserve the historical drafts. Public Markdown manuscripts have been edited to remove internal notes and placeholders, repair conversion artifacts and align figure numbers with the final wiki. They retain historical scientific claims that are not necessarily adopted by the portfolio narrative. The [published team wiki](https://2025.igem.wiki/sjtu-biox-shanghai/model) is the primary presentation of the final results.

| Source issue | Treatment in the portfolio |
| --- | --- |
| `Overview.docx` says scFv 42 mM; `MODEL Part1.docx` says approximately 0.42 mM | Quantitative molecular secretion claim omitted pending verification of the model output and concentration units |
| Overview describes ODE work as “in vitro” | Curated description identifies it as computational/in-silico modeling |
| scFv and VHH names intermixed | Molecular scFv analysis and proposed VHH design are distinguished |
| Simulation described as confirming therapeutic efficacy | Curated text describes model predictions and avoids clinical efficacy claims |
| Recommendation of antibody injections in results draft | Omitted from curated narrative; the project does not provide patient treatment guidance |
| `Recurrence/config/settings.xml` references the Treatment rule-file directory | Portable execution uses the Recurrence project's own supplied `cell_rules.csv`; historical XML remains available in the linked GitLab source |
| Absolute paths from the original machines | Removed from published configuration and metadata copies; transformations are recorded in the source manifest |
| 60–62 OpenMP threads in historical settings | Portable runner defaults to 4 and accepts an explicit thread count |
| Active Studio file named `settings` or `settings.xml` while C++ defaults to `PhysiCell_settings.xml` | Runner explicitly selects and passes the active settings |
| Recurrence settings specify 2,880 minutes; archived sweep series end at 1,440 minutes | Documented as separate full-configuration and archived-data horizons |
| Wiki references movies and full substrate fields not included in GitLab | Embedded figures/animations and existing population data are included locally; the manuscripts now link additional animations and reference PDFs hosted on the official wiki |
| Draft swaps explanations of Studio `-c` and `-e` | Run guide gives `-c` for configuration and `-e` for executable |
| Original Word field codes inside equations and plain-text LaTeX | RNAII/RNAIII field artifacts and multiline table cells repaired for GitHub; original DOCX retained locally |
| AIP-C inhibition half-max reported as 0.667 in the wiki; active Recurrence, Treatment and Recurrence_3D rule CSVs use 0.001 | Retained as a source discrepancy pending original calibration records; no biological parameters changed during portfolio editing |
| Recurrence_3D description refers to direct VHH suppression of cytokines; no direct VHH-to-cytokine rule was located in its supplied active rules or custom reaction implementation | Full historical 3D mechanism remains unverified from the supplied source; indirect toxin-mediated effects are possible |

The reference lists in historical drafts are preserved as supplied. This packaging work does not constitute a bibliographic audit of every biological parameter source. The provenance manifest and upstream license notices allow readers to distinguish original materials from the 2026 portfolio additions.

The October 2026 rendering check recovered three lost Word crop settings. The public manuscripts now reference unmodified published wiki images for Figures 3, 4 and 15, while retaining the original embedded media. Raw HTML image dimensions were replaced with Markdown image references; standalone LaTeX equations were given math blocks.
