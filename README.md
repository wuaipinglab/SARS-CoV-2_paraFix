# Test parallel and fixation sites

To test how much the fixation and parallel mutations proposed in evolution theory reflects the future and their predictivity power if there is any.

## Dependencies to run this project

Python packages: [`biopython`](https://biopython.org/wiki/Download), [`pandas`](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html), [`matplotlib`](https://matplotlib.org/stable/users/installing.html). Here is one of the ways to install them.

```bash
pip install biopython pandas matplotlib
```

R packages: [`phangorn`](https://cran.r-project.org/package=phangorn), [`VennDiagram`](https://CRAN.R-project.org/package=VennDiagram), [`gridExtra`](https://CRAN.R-project.org/package=gridExtra), [`svglite`](https://CRAN.R-project.org/package=svglite) and [`jsonlite`](https://cran.r-project.org/package=jsonlite).

```r
install.packages(c("phangorn", "VennDiagram", "gridExtra", "svglite", "jsonlite"))
```

R Bioconductor package: [`sitePath`](https://bioconductor.org/packages/sitePath/) and [`trackViewer`](	https://bioconductor.org/packages/trackViewer/). Here is one of the ways to install it.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("sitePath", "trackViewer"))
```

External packages: [`hyphy`](https://github.com/veg/hyphy) and [`muscle`](https://www.drive5.com/muscle/). If you use [Bioconda](https://bioconda.github.io/), the installation can be easily done by
```bash
conda install hyphy muscle
```
