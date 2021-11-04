# Test parallel and fixation sites

To test how much the fixation and parallel mutations proposed in evolution theory reflects the future and their predictivity power if there is any.

## Dependencies to run this project

Python packages: [`biopython`](https://biopython.org/wiki/Download), [`pandas`](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html), [`matplotlib`](https://matplotlib.org/stable/users/installing.html). Here is one of the ways to install them.

```bash
pip install biopython pandas matplotlib
```

R packages: [`phangorn`](https://cran.r-project.org/package=phangorn), [`jsonlite`](https://cran.r-project.org/package=jsonlite).

```r
install.packages(c("phangorn", "jsonlite"))
```

R Bioconductor package: [`sitePath`](https://bioconductor.org/packages/sitePath/). Here is one of the ways to install it.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sitePath")
```

External packages: [`hyphy`](https://github.com/veg/hyphy). If you use [Bioconda](https://bioconda.github.io/), the installation can be easily done by
```bash
conda install hyphy
```
