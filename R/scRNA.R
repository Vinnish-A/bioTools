
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
#' cluster2major(cluster_A = c(0, 1), cluster_B = c(2, 3), cluster_C = NULL, nCluster_ = 5)
#'
#' @export
cluster2major = function(..., nCluster_, remain2_ = NULL) {

  exprs_ = exprs(...)
  exprs_ = lapply(exprs_, eval)

  clusters_ = 1:nCluster_ - 1
  existing_ = unlist(exprs_)

  assertthat::assert_that(all(existing_ %in% clusters_))

  null_ = setdiff(clusters_, unlist(exprs_))
  if (!all(clusters_ %in% existing_)) {

    if (is.null(remain2_)) {
      exprs_[[which(sapply(exprs_, is.null))]] = null_
    } else {
      exprs_[[remain2_]] = c(exprs_[[remain2_]], null_)
    }

  }

  cluster_ = as.character(unlist(exprs_))
  cell_ = factor(rep(names(exprs_), sapply(exprs_, length)), levels = names(exprs_))

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
#' cl2ma = cluster2major(cluster_A = c(0, 1), cluster_B = c(2, 3), cluster_C = NULL, nCluster_ = 5)
#' major2minor(cluster_A1 = 0, cluster_A2 = 1, cl2ma_ = cl2ma, attention_ = "cluster_A")
#'
#' @export
major2minor = function(..., cl2ma_, attention_, remain2_ = NULL) {

  exprs_ = exprs(...)
  exprs_ = lapply(exprs_, eval)

  remain_ = as.numeric(names(cl2ma_[cl2ma_ == attention_]))
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
  cell_ = factor(rep(names(exprs_), sapply(exprs_, length)), levels = names(exprs_))

  table_min_ = setNames(cell_, cluster_)

  levels_ = as.list(levels(cl2ma_))
  levels_[[which(sapply(levels_, \(level__) level__ == attention_))]] = levels(cell_)
  levels_ = unlist(levels_)

  table_ = factor(c(cl2ma_[cl2ma_ != attention_], table_min_), levels = levels_)

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
#'
#' @export
mergeSeurat = function(..., lst_ = NULL) {

  if (is.null(lst_)) lst_ = list(...)

  x_ = lst_[[1]]
  y_ = lst_[-1]

  merge(x_, y_, add.cell.ids = T)

}

# unstable
get_expr = function(seu_, fea_, rename2expr_ = F) {

  if ('cell' %in% colnames(seu_@meta.data)) {

    res_ =  inner_join(
      as_tibble(FetchData(seu_, fea_), rownames = 'cell'),
      as_tibble(seu_@meta.data)
    )

  } else {

    res_ =  inner_join(
      as_tibble(FetchData(seu_, fea_), rownames = 'cell'),
      as_tibble(seu_@meta.data, rownames = 'cell')
    )

  }

  if (length(fea_) == 1 & rename2expr_) {
    colnames(res_)[[2]] = 'expression'
  }

  return(res_)

}

# unstable
get_umap = function(seu_) {

  if ('cell' %in% colnames(seu_@meta.data)) {

    res_ =  inner_join(
      as_tibble(seu_@reductions$umap@cell.embeddings, rownames = 'cell'),
      as_tibble(seu_@meta.data)
    )

  } else {

    res_ =  inner_join(
      as_tibble(seu_@reductions$umap@cell.embeddings, rownames = 'cell'),
      as_tibble(seu_@meta.data, rownames = 'cell')
    )

  }

  return(res_)

}

# unstable
drop_outlier_each = function(dp_, x_ = 'umap_1', y_ = 'umap_2', distance_ = 0.5, nNeighbor_ = 5) {

  checkReliance('dbscan')

  cols_ = c(x_, y_)

  vec_group_ = dp_ |>
    select(all_of(cols_)) |>
    dbscan::dbscan(eps = distance_, minPts = nNeighbor_) |>
    _[['cluster']]

  dp_ |>
    filter(vec_group_ == 1)

}

# unstable
drop_outlier = function(dp_, group_, x_ = 'umap_1', y_ = 'umap_2', distance_ = 0.5, nNeighbor_ = 5) {

  lst_dp_ = split(dp_, dp_[[group_]])

  lst_res_ = map(lst_dp_, drop_outlier_each, x_ = x_, y_ = y_, distance_ = distance_, nNeighbor_ = nNeighbor_)

  bind_rows(lst_res_)

}


