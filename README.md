# MoussaMap

MoussaMap is an R package for building interactive heatmaps from an expression matrix and cell metadata.

## Install from GitHub

```r
install.packages("pak")
pak::pkg_install("Fnorgh/MoussaMap")
```

Or with `remotes`:

```r
install.packages("remotes")
remotes::install_github("Fnorgh/MoussaMap")
```

## Use the package

```r
library(MoussaMap)

result <- MoussaMap(
	expression_matrix = "expr_matrix.csv",
	metadata = "metadata.csv",
	normalization = "none"
)

result$plot
```
