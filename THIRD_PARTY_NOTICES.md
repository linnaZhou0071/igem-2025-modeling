# Third-party software and attribution

## PhysiCell 1.14.2 and BioFVM

Source components are copied from the local PhysiCell 1.14.2 installation into `vendor/PhysiCell`. Original source copyright/license headers are retained. Component licenses are preserved in [vendor/PhysiCell/licenses](vendor/PhysiCell/licenses), including PhysiCell, BioFVM and pugixml. The original project-specific C++ files also retain their upstream notices.

Primary engine citations:

- Ghaffarizadeh A, Heiland R, Friedman SH, Mumenthaler SM, Macklin P. **PhysiCell: An open source physics-based cell simulator for 3-D multicellular systems.** *PLoS Computational Biology* 14(2), e1005991 (2018). [DOI](https://doi.org/10.1371/journal.pcbi.1005991).
- Ghaffarizadeh A, Friedman SH, Macklin P. **BioFVM: An efficient, parallelized diffusive transport solver for 3-D biological simulations.** *Bioinformatics* 32(8), 1256–1258 (2016). [DOI](https://doi.org/10.1093/bioinformatics/btv730).

Additional citation information is retained in [CITATION.txt](vendor/PhysiCell/CITATION.txt) and the original recurrence [ALL_CITATIONS.txt](models/Recurrence/ALL_CITATIONS.txt).

## Other tools

The documented molecular model used MATLAB SimBiology, which is separately licensed proprietary software and is not included here. The documentation discusses use of PhysiCell Studio, whose GUI is not bundled. Optional Python analysis dependencies are listed in `requirements.txt` and retain their own licenses.

No new repository-wide license is assigned to the original iGEM team materials, wiki figures or helper scripts. Existing notices continue to apply. Displaying an author's contribution does not transfer ownership of upstream or collaborative work.
