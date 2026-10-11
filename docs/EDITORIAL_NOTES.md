# Editorial and evidence notes

The main README and research overview are curated portfolio descriptions. The Word files and their Markdown conversions preserve historical draft content, including editing notes and stronger claims that are not adopted by the portfolio narrative.

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
| Wiki references movies and full substrate fields not included in GitLab | Only embedded figures/animations and existing population data are presented as available |
| Draft swaps explanations of Studio `-c` and `-e` | Run guide gives `-c` for configuration and `-e` for executable |
| Original Word field codes inside equations and plain-text LaTeX | Converted files may preserve equation artifacts; original DOCX is retained for reference |

The reference lists in historical drafts are preserved as supplied. This packaging work does not constitute a bibliographic audit of every biological parameter source. The provenance manifest and upstream license notices allow readers to distinguish original materials from the 2026 portfolio additions.
