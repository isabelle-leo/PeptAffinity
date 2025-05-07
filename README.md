# **PeptOlink**

*A Shiny app for exploring the relationship between MS–derived peptide intensities and paired Olink assays.*

Read the preprint here: https://www.researchsquare.com/article/rs-6501601/v1

PeptOlink maps peptides detected in plasma proteomics data by MS onto their parent protein sequences, annotates isoforms and features (InterPro/Prosite domains, structural features), and visualizes correlations between peptide abundances and Olink assay measurements. This helps identify regions of a protein driving Olink assay signals, and might highlight structural or domain hotspots of functional relevance in the plasma proteome.

**Getting started**

The app can be accessed at: https://peptolink.serve.scilifelab.se/

**Using the app code**

The app is hosted by scilifelab serve. All datasets used in the app are available in the repository data folder. If you want to run/adapt the app yourself, you can clone the repository from the latest stable tag version, set up dependencies in the renv.lock file, and deploy a docker using the dockerfile. You should supply your own docker and google analytics credentials. A data preprocessing script is also provided in the preprocessing/ folder, this can be used e.g. to update annotations with a new database version. Databases used in the app were last obtained: April 16, 2025.

The app is organized by: setup section, functions section, app section. Constants and style preferences are defined in the first section, all major and plotting functions are defined in the second section, and the ui and server of the shiny app are established in the third.

**User guide**

Filters panel:

* Correlation Threshold (≤): Upper limit for peptide‑Olink correlation

* Spread Threshold (≥): Minimum variability (SD/IQR) across samples

* Peptides (≥): Minimum number of distinct peptides per gene ID

* Isoforms (≥): Minimum number of protein isoforms detected

* Click “Clear All” to reset filters to default values
* The filter panel updates values automatically
* The tab retracts when you click again

Feature tab:

* Select Gene & Isoform: Choose a gene symbol and its UniProt isoform

* Domain Source: Toggle between InterPro and Prosite annotations

* Interactive Heatmap: Rows show peptides mapped along sequence; color scale indicates correlation value between MS abundance and Olink

* Tooltips: Hover on any tile to see residue position, correlation, and feature annotation

Structural Tab: *coming soon*

Summary tab:

* Use this tab to guide your selections in the filter panel and to see how much of the dataset meets your parameters

* Violin & Jitter: Distribution of all individual peptide correlations across the filtered set of gene IDs.

* Scatter plots: Mean vs. SD and median vs. IQR plots for gene ID‑level statistics

# Our preprint (please cite if you use our app): [https://www.researchsquare.com/article/rs-6501601/v1](https://doi.org/10.21203/rs.3.rs-6501601/v1)

© 2025 Isabelle Leo Noora Sissala Haris Babačić
Apache 2.0 license


