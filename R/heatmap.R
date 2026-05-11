#' Select variable features after averaging
#'
#' Internal helper used by `create_heatmap()` when `top_var_features` is set.
#'
#' @keywords internal
select_variable_features_after_averaging <- function(expr_matrix,
                                                    nfeatures = 2000,
                                                    selection.method = "vst") {
  expr_matrix <- as.matrix(expr_matrix)

  seu <- Seurat::CreateSeuratObject(counts = expr_matrix)
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(
    seu,
    selection.method = selection.method,
    nfeatures = min(nfeatures, nrow(expr_matrix)),
    verbose = FALSE
  )

  Seurat::VariableFeatures(seu)
}

read_optional_csv_vector <- function(path, header = FALSE) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }

  values <- read.csv(path, header = header, stringsAsFactors = FALSE)[[1]]
  values <- unique(trimws(as.character(values)))
  values <- values[values != ""]

  if (length(values) == 0) {
    return(NULL)
  }

  values
}

read_optional_group_level <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }

  group_info <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)
  chosen_level <- trimws(as.character(group_info[[1]][1]))

  if (is.na(chosen_level) || chosen_level == "") {
    return(NULL)
  }

  chosen_level
}

#' Create an interactive heatmap from expression and metadata CSV files
#'
#' @param expression_matrix Path to the expression CSV file.
#' @param metadata Path to the metadata CSV file.
#' @param gene_selection Optional path to a single-column CSV of genes to keep.
#' @param cell_selection Optional path to a single-column CSV of cells to keep.
#' @param group_by Optional path to a CSV containing the finest grouping level.
#' @param group_order Optional path to a CSV specifying metadata column order.
#' @param group_by_level Optional finest-level metadata column name.
#' @param top_var_features Optional number of variable features to keep.
#' @param varfeat_method Variable feature selection method passed to Seurat.
#' @param colorscale Heatmap color scale.
#' @param normalization Normalization mode.
#' @param width Plot width.
#' @param height Plot height.
#'
#' @return A list with `plot` and, on failure, `error`.
#' @export
create_heatmap <- function(expression_matrix,
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
                          height = 800) {
  tryCatch({
    if (!file.exists(expression_matrix)) {
      stop("Expression matrix CSV not found: ", expression_matrix)
    }
    if (!file.exists(metadata)) {
      stop("Metadata CSV not found: ", metadata)
    }

    expr_temp <- read.csv(expression_matrix, stringsAsFactors = FALSE, check.names = FALSE)
    gene_names <- make.unique(trimws(as.character(expr_temp[, 1])), sep = "_")
    rownames(expr_temp) <- gene_names
    expression_data <- expr_temp[, -1, drop = FALSE]
    metadata_data <- read.csv(metadata, stringsAsFactors = FALSE, check.names = FALSE)

    if (is.null(names(metadata_data)[1]) || names(metadata_data)[1] == "") {
      names(metadata_data)[1] <- "cell_id"
    }
    id_col <- names(metadata_data)[1]
    metadata_cols <- names(metadata_data)[-1]

    if (length(metadata_cols) == 0) {
      stop("Metadata must contain at least one grouping column after the cell ID column.")
    }

    common_cells <- intersect(colnames(expression_data), metadata_data[[id_col]])
    if (length(common_cells) == 0) {
      stop("No matching cells found between expression matrix and metadata.")
    }

    cell_list <- read_optional_csv_vector(cell_selection, header = FALSE)
    if (!is.null(cell_list)) {
      common_cells <- intersect(common_cells, cell_list)
    }

    expression_data <- expression_data[, !duplicated(colnames(expression_data)), drop = FALSE]
    metadata_data <- metadata_data[!duplicated(metadata_data[[id_col]]), , drop = FALSE]
    expression_data <- expression_data[, common_cells, drop = FALSE]
    metadata_data <- metadata_data[metadata_data[[id_col]] %in% common_cells, , drop = FALSE]
    idx <- match(common_cells, metadata_data[[id_col]])
    metadata_data <- metadata_data[idx, , drop = FALSE]

    gene_list <- read_optional_csv_vector(gene_selection, header = FALSE)
    if (!is.null(gene_list)) {
      expression_data <- expression_data[rownames(expression_data) %in% gene_list, , drop = FALSE]
    }

    if (normalization == "log2" || normalization == "log2_scale") {
      expression_data <- log2(expression_data + 1)
    }

    chosen_level <- read_optional_group_level(group_by)
    if (!is.null(chosen_level)) {
      group_by_level <- chosen_level
    }

    if (!is.null(group_order) && file.exists(group_order)) {
      order_info <- read.csv(group_order, header = TRUE, stringsAsFactors = FALSE)
      group_order <- unique(na.omit(trimws(unlist(order_info))))
      group_order <- group_order[group_order != ""]

      if (length(group_order) > 0) {
        if ("class" %in% group_order && "class" %in% names(metadata_data)) {
          group_by_level <- "class"
        }
        metadata_cols <- intersect(group_order, metadata_cols)
      } else {
        warning("group_order was empty or invalid.")
        group_order <- NULL
      }
    }

    metadata_cols <- trimws(metadata_cols)

    if (!is.null(group_by_level) && group_by_level %in% metadata_cols) {
      finest_level <- group_by_level
    } else {
      finest_level <- tail(metadata_cols, 1)
    }

    finest_index <- match(finest_level, metadata_cols)
    if (is.na(finest_index) || finest_index < 1) {
      stop("Could not determine a valid grouping level from metadata columns.")
    }
    grouping_cols <- metadata_cols[seq_len(finest_index)]

    metadata_data$full_path <- apply(
      metadata_data[, grouping_cols, drop = FALSE],
      1,
      function(x) paste(x, collapse = "_")
    )

    cell_assignments <- setNames(metadata_data$full_path, metadata_data[[id_col]])

    combined_data <- expression_data |>
      tibble::rownames_to_column("gene") |>
      tidyr::pivot_longer(-gene, names_to = "cell", values_to = "expression") |>
      dplyr::mutate(cell_group = cell_assignments[cell]) |>
      dplyr::filter(!is.na(cell_group))

    cell_type_averages <- combined_data |>
      dplyr::group_by(gene, cell_group) |>
      dplyr::summarise(avg_expression = mean(expression, na.rm = TRUE), .groups = "drop")

    expr_matrix <- cell_type_averages |>
      tidyr::pivot_wider(names_from = cell_group, values_from = avg_expression, values_fill = 0) |>
      tibble::column_to_rownames("gene") |>
      as.matrix()

    if (any(is.na(expr_matrix)) || any(is.infinite(expr_matrix))) {
      expr_matrix[is.na(expr_matrix)] <- 0
      expr_matrix[is.infinite(expr_matrix)] <- 0
    }

    if (!is.null(top_var_features) && is.numeric(top_var_features) && top_var_features > 0) {
      var_genes <- select_variable_features_after_averaging(
        expr_matrix = expr_matrix,
        nfeatures = top_var_features,
        selection.method = varfeat_method
      )
      var_genes <- intersect(var_genes, rownames(expr_matrix))
      if (length(var_genes) == 0) {
        warning("top_var_features requested, but no VariableFeatures matched rownames(expr_matrix). Skipping.")
      } else {
        expr_matrix <- expr_matrix[var_genes, , drop = FALSE]
      }
    }

    if (normalization == "scale" || normalization == "log2_scale") {
      expr_matrix <- t(scale(t(expr_matrix)))
    } else if (normalization == "max") {
      max_vals <- apply(expr_matrix, 1, max, na.rm = TRUE)
      max_vals[max_vals == 0] <- 1
      expr_matrix <- expr_matrix / max_vals
    }

    col_groups <- colnames(expr_matrix)
    hierarchy_df <- data.frame(full_path = col_groups, stringsAsFactors = FALSE)

    for (i in seq_along(grouping_cols)) {
      hierarchy_df[[grouping_cols[i]]] <- sapply(strsplit(col_groups, "_"), function(x) {
        if (length(x) >= i) x[i] else NA
      })
    }

    col_side_colors <- list()
    legend_data <- list()

    hue_families <- list(
      c("#08306b", "#6baed6", "#2171b5", "#c6dbef", "#08519c", "#9ecae1", "#4292c6", "#deebf7"),
      c("#67000d", "#fb6a4a", "#a50f15", "#fcbba1", "#cb181d", "#fc9272", "#ef3b2c", "#fee0d2"),
      c("#00441b", "#74c476", "#238b45", "#c7e9c0", "#006d2c", "#a1d99b", "#41ab5d", "#e5f5e0"),
      c("#3f007d", "#9e9ac8", "#54278f", "#dadaeb", "#6a51a3", "#bcbddc", "#807dba", "#f2f0f7"),
      c("#7f2704", "#fd8d3c", "#a63603", "#fdd0a2", "#d94801", "#fdae6b", "#f16913", "#fee6ce"),
      c("#003c30", "#80cdc1", "#01665e", "#c7eae5", "#35978f", "#b2e2e2", "#5ab4ac", "#edf8fb")
    )

    for (level_idx in seq_along(grouping_cols)) {
      level_name <- grouping_cols[level_idx]
      level_values <- hierarchy_df[[level_name]]
      unique_vals <- unique(level_values[!is.na(level_values)])
      n_unique <- length(unique_vals)

      pal <- hue_families[[(level_idx - 1) %% length(hue_families) + 1]]
      colors <- if (n_unique <= length(pal)) {
        pal[seq_len(n_unique)]
      } else {
        grDevices::colorRampPalette(pal)(n_unique)
      }

      color_map <- setNames(colors, unique_vals)
      color_vals <- paste(level_values, color_map[level_values], sep = "||")
      col_side_colors[[level_name]] <- color_vals
      legend_data[[level_name]] <- data.frame(label = unique_vals, color = colors)
    }

    col_side_colors <- as.data.frame(col_side_colors, row.names = col_groups)
    col_side_colors <- col_side_colors[, rev(colnames(col_side_colors)), drop = FALSE]

    color_func <- switch(
      colorscale,
      viridis = viridis::viridis,
      plasma = viridis::plasma,
      Blues = function(n) RColorBrewer::brewer.pal(n, "Blues"),
      Reds = function(n) RColorBrewer::brewer.pal(n, "Reds"),
      viridis::viridis
    )

    if (!is.null(group_order)) {
      dendro_mode <- "row"
      show_dend <- c(TRUE, FALSE)
    } else {
      dendro_mode <- "both"
      show_dend <- c(TRUE, TRUE)
    }

    norm_label <- switch(
      normalization,
      "log2" = " (log2 normalized)",
      "scale" = " (scaled)",
      "log2_scale" = " (log2 + scaled)",
      "max" = " (max scaled)",
      "none" = " (no normalization)",
      ""
    )

    col_side_colors[] <- lapply(col_side_colors, function(x) as.character(x))
    heatmap_plot <- heatmaply::heatmaply(
      expr_matrix,
      colors = color_func(256),
      dendrogram = dendro_mode,
      show_dendrogram = show_dend,
      col_side_colors = col_side_colors,
      main = paste0("Expression Heatmap - grouped by ", finest_level, norm_label),
      xlab = "Groups",
      ylab = "Genes",
      fontsize_row = 8,
      fontsize_col = 6,
      hclust_method = "ward.D2",
      dist_method = "euclidean",
      width = width,
      height = height
    )

    annotation_order <- colnames(col_side_colors)
    plt <- heatmap_plot

    legend_entries <- list()
    seen <- character()

    for (level in annotation_order) {
      for (tr in plt$x$data) {
        if (!is.null(tr$fillcolor) && !is.null(tr$text)) {
          variable <- sub(".*variable: ([^<]+).*", "\\1", tr$text)
          value <- sub(".*value: ([^<]+).*", "\\1", tr$text)

          if (variable != level) {
            next
          }

          key <- paste(variable, value, tr$fillcolor)

          if (!(key %in% seen)) {
            seen <- c(seen, key)
            legend_entries[[length(legend_entries) + 1]] <- list(
              type = "scatter",
              mode = "markers",
              x = NA,
              y = NA,
              marker = list(size = 10, color = tr$fillcolor),
              name = paste(variable, "=", value),
              legendgroup = variable,
              showlegend = TRUE,
              hoverinfo = "skip"
            )
          }
        }
      }
    }

    for (i in seq_along(plt$x$data)) {
      tr <- plt$x$data[[i]]
      if (!is.null(tr$name) && grepl("^#[0-9A-Fa-f]{6}$", tr$name)) {
        plt$x$data[[i]]$showlegend <- FALSE
      }
    }

    plt$x$data <- c(plt$x$data, legend_entries)
    heatmap_plot <- plt

    new_data <- list()
    seen_legend <- character()

    for (tr in plt$x$data) {
      if (!is.null(tr$fillcolor) && !is.null(tr$text) && !is.null(tr$x)) {
        level <- sub(".*variable: ([^<]+).*", "\\1", tr$text)
        if (!level %in% grouping_cols) {
          new_data[[length(new_data) + 1]] <- tr
          next
        }

        group_name <- tr$x[1]
        idx <- match(group_name, hierarchy_df$full_path)
        if (is.na(idx)) {
          new_data[[length(new_data) + 1]] <- tr
          next
        }

        meta_value <- hierarchy_df[[level]][idx]
        tr$text <- paste0("variable: ", level, "<br>", "value: ", meta_value)

        legend_key <- paste(level, meta_value)
        tr$name <- meta_value
        tr$legendgroup <- level

        if (!(legend_key %in% seen_legend)) {
          tr$showlegend <- TRUE
          seen_legend <- c(seen_legend, legend_key)
        } else {
          tr$showlegend <- FALSE
        }

        new_data[[length(new_data) + 1]] <- tr
      } else {
        new_data[[length(new_data) + 1]] <- tr
      }
    }

    plt$x$data <- new_data
    heatmap_plot <- plt

    plt <- heatmap_plot
    for (i in seq_along(plt$x$data)) {
      tr <- plt$x$data[[i]]
      if (isTRUE(tr$showlegend) && !is.null(tr$name)) {
        if (grepl("\\|\\|#", tr$name)) {
          cleaned <- sub(".*=\\s*", "", tr$name)
          cleaned <- sub("\\|\\|#.*$", "", cleaned)
          level <- tr$legendgroup
          plt$x$data[[i]]$name <- paste0(level, ": ", cleaned)
        }
      }
    }

    heatmap_plot <- plt
    heatmap_plot$x$layout$legend <- list(
      orientation = "v",
      x = 1.22,
      y = 1,
      xanchor = "left",
      yanchor = "top",
      font = list(size = 12)
    )

    heatmap_plot$x$layout$margin <- list(
      l = 80,
      r = 260,
      b = 80,
      t = 80
    )

    heatmap_plot$x$layout$showlegend <- TRUE
    heatmap_plot$x$data <- c(heatmap_plot$x$data, list())

    plt <- heatmap_plot
    for (i in seq_along(plt$x$data)) {
      tr <- plt$x$data[[i]]
      if (!is.null(tr$x) && !all(is.na(tr$x))) {
        plt$x$data[[i]]$showlegend <- FALSE
      }
    }

    heatmap_plot <- plt
    list(plot = heatmap_plot)
  }, error = function(e) {
    traceback()
    list(plot = NULL, error = paste("Error:", e$message))
  })
}


