# **PeptAffinity**

*A Shiny app for exploring the relationship between MS–derived peptide intensities and paired Olink assays.*

Read the publication here: [DOI: 10.1038/s42004-025-01753-2](https://doi.org/10.1038/s42004-025-01753-2) and please cite if you use the app!

PeptAffinity maps peptides detected in plasma proteomics data by MS onto their parent protein sequences, annotates isoforms and features (InterPro/Prosite domains, structural features), and visualizes correlations between peptide abundances and Olink assay measurements. This helps identifying regions of a protein that correlate to Olink assay measurements, and might highlight structural or domain hotspots of functional relevance in the plasma proteome.

**Getting started**

The app can be accessed at: https://peptaffinity.serve.scilifelab.se/

**Using the app code**

The app is hosted by scilifelab serve. All datasets used in the app are available in the repository data folder. If you want to run/adapt the app yourself, you can clone the repository from the latest stable tag version, set up dependencies in the renv.lock file, and deploy a docker using the dockerfile. You should supply your own docker and google analytics credentials. A data preprocessing script is also provided in the preprocessing/ folder, this can be used e.g. to update protein annotations with a new database version. Proteome databases used in the app were last obtained: April 16, 2025. The structural models used in the app were updated: October 20, 2025.

The app is organized by: setup section, functions section, app section. Constants and style preferences are defined in the first section, all major and plotting functions are defined in the second section, and the ui and server of the shiny app are established in the third.

**User guide**

Protein Filters panel:

* Peptide Correlation: Range of peptide‑Olink correlations (based on the mean, median, or center values of peptide correlation per protein)

* Peptide Correlation Spread: Range for the variance across samples (based on standard deviation - SD, interquartile range - IQR, or range of peptide correlation per protein)

* Number of Peptides (≥): Minimum number of distinct peptides per gene ID

* Number of Isoforms (≥): Minimum number of protein isoforms detected

* Click “Clear All” to reset filters to include all valid peptide correlations
* The filter panel updates values automatically
* The tab retracts when you click again

Peptide Filters panel:

* Samples with peptide: Define sample detection criteria for peptides plotted (does not impact gene symbol filtering or summary statistics)

Home tab:

* Here you can find the tabs' structure for the app

Sequence tab:

* Select Gene & Isoform: Choose a gene symbol and its UniProt isoform

* Domain Source: Toggle between InterPro and Prosite annotations

* Interactive Heatmap of the protein: Rows show peptides mapped along sequence; color scale indicates correlation between MS abundance and Olink NPX values

* Tooltips: Hover on any tile to see residue position, correlation, and feature annotation

Structure tab:

* Click "Spin" in order to rotate the structure 

* Interactive AlphaFold structure of the protein: the color scale indicates median peptide correlation between MS abundance and Olink NPX values

All proteins tab:

* Use this tab to guide your selections in the filter panel and to see how much of the dataset meets your parameters

* Scatter plots: Mean vs. SD and center vs. range plots for gene ID‑level statistics

# Our paper (please cite if you use our app): [Sissala, N., Babačić, H., Leo, I.R. et al. Comparative evaluation of Olink Explore 3072 and mass spectrometry with peptide fractionation for plasma proteomics. Commun Chem 8, 327 (2025).](https://doi.org/10.1038/s42004-025-01753-2)
Our previous preprint: [https://www.researchsquare.com/article/rs-6501601/v1](https://doi.org/10.21203/rs.3.rs-6501601/v1)

© 2025 Isabelle Leo Noora Sissala Haris Babačić
Apache 2.0 license


