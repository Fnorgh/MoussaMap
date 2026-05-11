# MoussaMap

MoussaMap is an R package for building interactive heatmaps from an expression matrix and cell metadata.

It is designed to:

- read a gene-by-cell expression matrix and matching cell metadata
- optionally filter genes and cells from CSV selection files
- group cells into hierarchical paths such as `class -> subclass -> cluster`
- average expression within each group
- optionally normalize or scale the grouped matrix
- return an interactive `heatmaply` plot with annotation bars and a cleaned legend

## Install from GitHub

```r
install.packages("pak")
pak::pkg_install("Fnorgh/MoussaMap")
```

Or with `remotes`:

```r
install.packages("remotes")
remotes::install_github("Fnorgh/MoussaMap", upgrade = "never")
```


## Quick Start

```r
library(MoussaMap)

result <- MoussaMap(
  expression_matrix = "expr_matrix.csv",
  metadata = "metadata.csv",
  normalization = "none"
)

result$plot
```

The function returns a list with:

- `result$plot`: the interactive heatmap
- `result$error`: an error message if something failed

## Required Input Files

You must provide two CSV files.

### 1. Expression Matrix CSV

Passed to `expression_matrix`.

Requirements:

- rows are genes
- columns are cells
- the first column contains gene names

Example:

```csv
gene,cellA,cellB,cellC
MS4A1,0,1.2,0.3
LST1,4.5,0.1,2.2
IL7R,2.0,0.0,1.1
```

### 2. Metadata CSV

Passed to `metadata`.

Requirements:

- the first column contains cell IDs
- those cell IDs must match the expression matrix column names
- the remaining columns are hierarchy levels ordered from broad to fine

Example:

```csv
cell_id,class,subclass,cluster
cellA,Immune,Tcell,Tcell_1
cellB,Immune,Bcell,Bcell_2
cellC,Glia,Microglia,Micro_3
```

Important details:

- the first metadata column is treated as the cell ID column even if it has a different name
- if the first metadata column has no name, the function renames it to `cell_id`
- all metadata columns after the first are treated as grouping levels

## Optional Input Files

### `gene_selection`

Optional path to a one-column CSV listing genes to keep.

Example:

```csv
MS4A1
LST1
IL7R
```

### `cell_selection`

Optional path to a one-column CSV listing cells to keep.

Example:

```csv
cellA
cellC
```

### `group_by`

Optional path to a CSV containing the finest grouping level to use.

Example:

```csv
level
cluster
```

### `group_order`

Optional path to a CSV listing the metadata columns to use and their order.

Example:

```csv
class
subclass
cluster
```

If `group_order` is provided, the heatmap keeps row clustering and hides the column dendrogram.

## Main Arguments

```r
MoussaMap(
  expression_matrix,
  metadata,
  gene_selection = NULL,
  cell_selection = NULL,
  group_by = NULL,
  group_order = NULL,
  group_by_level = NULL,
  top_var_features = NULL,
  varfeat_method = "vst",
  colorscale = "viridis",
  normalization = "none",
  width = 1000,
  height = 800
)
```

### Advanced options

- `top_var_features`: keep only the top variable genes after averaging
- `varfeat_method`: method passed to Seurat for variable feature selection
- `group_by_level`: directly set the finest grouping column without using a file

## Normalization Options

Use the `normalization` argument.

- `"none"`: no transform
- `"log2"`: apply `log2(expr + 1)` before averaging
- `"scale"`: z-score each gene after averaging
- `"log2_scale"`: log-transform before averaging, then z-score after averaging
- `"max"`: divide each gene row by its maximum value after averaging

Notes:

- scaling is done row-wise by gene
- `NA` and `Inf` values are replaced with `0` before plotting

## Color Options

Use the `colorscale` argument.

- `"viridis"` default
- `"plasma"`
- `"Blues"`
- `"Reds"`

These control the heatmap intensity colors for expression values.

## What The Plot Shows

The output plot includes:

- clustered genes
- grouped cell columns
- annotation bars for each metadata hierarchy level used
- cleaned legend labels so users see category names instead of hex colors
- hover text showing annotation variable and value

## Example With More Options

```r
library(MoussaMap)

result <- MoussaMap(
  expression_matrix = "expr_matrix.csv",
  metadata = "metadata.csv",
  gene_selection = "selected_genes.csv",
  cell_selection = "selected_cells.csv",
  group_by = "group_by_level.csv",
  group_order = "group_order.csv",
  top_var_features = 200,
  varfeat_method = "vst",
  colorscale = "viridis",
  normalization = "log2_scale",
  width = 1200,
  height = 900
)

result$plot
```


### Or very basic

```r
library(MoussaMap)
result <- MoussaMap(
  expression_matrix = "expr_matrix.csv",
  metadata = "metadata.csv",
)
result$plot
```


