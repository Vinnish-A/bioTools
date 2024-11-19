
# scRNA -------------------------------------------------------------------


## Annotation ----

#' Map clusters to cell types
#'
#' Assigns cluster numbers to cell type labels. Unassigned clusters are grouped into the remaining group.
#'
#' @param ... Expressions that map cluster numbers to cell types.
#' @param nCluster_ Integer. The total number of clusters.
#' @param remain2_ Integer. The target group for unassigned clusters. Defaults to the first empty group if NULL.
#'
#' @return A named vector where keys are cluster numbers, and values are corresponding cell type labels.
#'
#' @examples
#' cluster2cell(cluster_A = c(0, 1), cluster_B = c(2, 3), nCluster_ = 5)
cluster2major = function(..., nCluster_, remain2_ = NULL) {

  exprs_ = exprs(...)
  exprs_ = lapply(exprs_, eval)

  clusters_ = 1:nCluster_ - 1
  existing_ = unlist(exprs_)

  assertthat::assert_that(all(existing_ %in% clusters_))

  null_ = setdiff(clusters_, unlist(exprs_))
  if (!all(sort(existing_) == clusters_)) {

    if (is.null(remain2_)) {
      exprs_[[which(sapply(exprs_, is.null))]] = null_
    } else {
      exprs_[[remain2_]] = c(exprs_[[remain2_]], null_)
    }

  }

  cluster_ = as.character(unlist(exprs_))
  cell_ = rep(names(exprs_), sapply(exprs_, length))

  table_ = setNames(cell_, cluster_)

  return(table_)

}

#' Update cell type assignments
#'
#' Updates the mapping of clusters to cell types based on the given focus cell type.
#' Unspecified clusters are added to the target group.
#'
#' @param ... Expressions that map cluster numbers to target groups.
#' @param cl2ma_ Named vector. Current mapping of clusters to cell types.
#' @param attention_ Character. The target cell type to focus on.
#' @param remain2_ Character. The target group for unassigned clusters. Defaults to be NULL.
#'
#' @return An updated named vector of cell type assignments.
#'
#' @importFrom assertthat assert_that
#'
#' @examples
#' cell2min(cluster_X = c(0, 1), cl2c_ = c("0" = "TypeA", "1" = "TypeA", "2" = "TypeB"), attention_ = "TypeA")
major2minor = function(..., cl2ma_, attention_, remain2_ = NULL) {

  exprs_ = exprs(...)
  exprs_ = lapply(exprs_, eval)

  remain_ = as.numeric(names(cl2c_[cl2c_ == attention_]))
  existing_ = unlist(exprs_)

  assertthat::assert_that(all(existing_ %in% remain_))

  null_ = setdiff(remain_, unlist(exprs_))

  if (!all(sort(existing_) == remain_)) {

    if (is.null(remain2_)) {
      exprs_[[which(sapply(exprs_, is.null))]] = null_
    } else {
      exprs_[[remain2_]] = c(exprs_[[remain2_]], null_)
    }

  }


  cluster_ = as.character(unlist(exprs_))
  cell_ = rep(names(exprs_), sapply(exprs_, length))

  table_min_ = setNames(cell_, cluster_)

  table_ = c(cl2c_[cl2c_ != attention_], table_min_)

  return(table_)

}

## manipulate ----

#' Merge Seurat objects
#'
#' Combines multiple Seurat objects into a single Seurat object.
#'
#' @param ... One or more Seurat objects to merge.
#' @param lst_ List. A list of Seurat objects to merge. Defaults to NULL, which uses the objects passed through `...`.
#'
#' @return A merged Seurat object.
#'
#' @examples
#' mergeSeurat(obj1, obj2, obj3)
mergeSeurat = function(..., lst_ = NULL) {

  if (is.null(lst_)) lst_ = list(...)

  x_ = lst_[[1]]
  y_ = lst_[-1]

  merge(x_, y_, add.cell.ids = T)

}

## simple_vis ----

#' Plot UMAP visualization
#'
#' Creates a UMAP plot with points colored by the specified grouping variable.
#'
#' @param df_ Data frame. Contains UMAP coordinates and grouping information.
#' @param x_ Character. The column name for the UMAP x-axis.
#' @param y_ Character. The column name for the UMAP y-axis.
#' @param colorBy_ Character. The column name for the grouping variable used for coloring points.
#'
#' @return A ggplot object displaying the UMAP visualization.
#'
#' @importFrom grid arrow
#'
#' @examples
#' plot_umap(df_, x_ = "UMAP_1", y_ = "UMAP_2", colorBy_ = "CellType")
plot_umap = function(df_, x_, y_, colorBy_) {

  x_ = sym(x_)
  y_ = sym(y_)
  colorBy_ = sym(colorBy_)

  axis_ = guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(4, "cm")
  )

  ggplot() +
    ggunchull::stat_unchull(aes(!!x_, !!y_, color = !!colorBy_, fill = !!colorBy_), alpha = 0.2, size = 1, lty = 2, qval = 0.8, delta = 1) +
    geom_point(aes(!!x_, !!y_, color = !!colorBy_), size = 1) +
    # geom_label(data = dataPlotUMAPText, aes(x, y, label = cluster, color = cluster)) +
    guides(fill = "none", x = axis_, y = axis_) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(color = "", title = "UMAP by CellType") +
    scale_fill_manual(values = color_macaron) +
    scale_color_manual(values = color_macaron) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.key = element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom",
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_text(hjust = 0.05, face = "italic"),
      axis.line.x = element_line(arrow = arrow(type = "open", length = unit(0.5, "cm"))),
      axis.line.y = element_line(arrow = arrow(type = "open", length = unit(0.5, "cm")))
    )

}


