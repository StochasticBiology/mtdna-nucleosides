# mtdna-nucleosides

Analysis of data for nucleoside treatment of control and POLG cell lines, and miscellaneous data files from and used in the experiments.

`analysis.R` reads the various datafiles and performs hypothesis testing, LMM analysis, and visualisation.

The dataset for the last part of the analysis is >25MB in original Excel format. `reduced-data.csv` is a condensed version of the original raw data, containing everything needed for this analysis, to be compatible with Github's file size rules. `full-data.csv.zip` is a zipped CSV-format version of the original dataset.

In `raw-data` there is a collection of Excel files corresponding to in-cell observations, mass spectrometry data, and details of primers used. The experimental setups and these files' connection to manuscript figures are described in `raw-data/Summary of files and runs.xlsx`.
