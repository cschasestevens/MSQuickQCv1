# MSQuickQC v1.0

Performs rapid quality control assessment of single or paired mass spectrometry datasets.

## Description

Untargeted metabolomics data must be technically robust to derive biological findings that are not masked by variance introduced by sample preparation and data acquisition. Determining and addressing potential sources of technical variance prior to statistical analysis is essential to perform unbiased data analysis and interpretation. MSQuickQC is a R Shiny app that screens mass spectrometry-based datasets to quantify sources of sample technical error and variation. The app utilizes a generic input format that enables analysis of both untargeted and targeted datasets, all major chromatography platforms and ionization modes, and multiple file types. Additionally, the app provides the option for side-by-side comparisons of datasets to visualize factors such as data normalization performance and instrument-to-instrument variability. The main output of the app is a reproducible report that includes a summary several commonly used metrics for assessing overall data quality.

## Getting Started

### Dependencies (if running on local machine)

* Windows 10-11 or macOS 12.7.1 or higher
* R version 4.3.1 or higher (https://cran.r-project.org/)
* RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* R-packages (downloaded from CRAN unless otherwise specified):
    * shinyBS
    * DT
    * shinyFeedback
    * shiny
    * IMIFA
    * ComplexHeatmap (Bioconductor)
    * readxl
    * ggplot2
    * purrr
    * ggpubr
    * ggrepel
    * RColorBrewer
    * circlize
    * corrplot
    * magrittr
    * viridis
    * patchwork
    * kableExtra

### Installing

* MSQuickQC is hosted under a free subscription at https://chasestevens.shinyapps.io/MSQuickQCv1/ which supports report generation for a **single file** due to memory limitations.
* Alternatively, the app can be run locally by downloading this repository or cloning using the following command in the console: ` git clone https://github.com/cschasestevens/MSQuickQCv1.git `
* For running local instances of MSQuickQC, ensure that R, RStudio, and all dependencies are installed. For users unfamiliar with R, packages can be downloaded by opening the included **install_packages.R** file in RStudio, highlighting all lines of the script, and pressing **Ctrl+Enter**.

### Executing program

* Once all dependencies have been installed, open **app.R** in RStudio and click the **Run App** button above the **top-left** pane.
* Follow the instructions included in the app UI to specify input data parameters and generate a report.

## Help

* All quality control names and identifier used for internal standards (if applicable) **must** match the names included in the input data. If these data are not included, simply type **NA** in each field.
* For reports including two datasets, ensure that entries are separated by a comma, without a space between entries (i.e. **Entry1,Entry2**). The first entry corresponds to input parameters of input file **1**, while the second entry corresponds to input parameters of input file **2**.

## Authors

Nathanial Chase Stevens, PhD
Nathanial_Stevens@med.unc.edu
cschasestevens@gmail.com
https://www.linkedin.com/in/nathanial-chase-stevens-phd-08775180/

## Version History

* 1.0
    * Initial Release

## License

This project is licensed under the GNU General Public License Version 3 - see the LICENSE.md file for details

## Acknowledgments

* IMIFA package: Murphy, K., Viroli, C., and Gormley, I. C. (2020). Infinite mixtures of infinite factor analysers. Bayesian Analysis, 15(3): 937–963. <doi:10.1214/19-BA1179>
* ComplexHeatmap package: Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics. <doi:10.1093/bioinformatics/btw313>
* Circlize package: Gu, Z. circlize implements and enhances circular visualization in R. Bioinformatics 2014.
