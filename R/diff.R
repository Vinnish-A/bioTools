
# diff --------------------------------------------------------------------


## check_input ----

checkGroup = function(group_) {

  assert_that(is.factor(group_))
  assert_that(levels(group_) == 2)

}

## DA ----

#' Differential Expression Analysis using Limma
#'
#' This function performs differential expression analysis using the Limma package.
#'
#' @param normal_matrix A numeric matrix of normalized expression values, where rows are genes and columns are samples.
#' @param group A factor vector indicating the group assignments for the samples.
#'
#' @return A tibble containing the differential expression results with columns:
#'   - `symbol`: Gene identifiers.
#'   - `logFC`: Log2 fold change.
#'   - `AveExpr`: Average expression level.
#'   - `t`: Moderated t-statistic.
#'   - `p.value`: Raw p-value.
#'   - `padj`: Adjusted p-value (FDR).
#'
#' @export
diff_limma = function(normal_matrix, group) {

  checkGroup(group)

  library(limma)

  design = model.matrix(~0 + group)
  colnames(design) = c('group1', 'group2')

  contrast = makeContrasts(group2 - group1, levels = design)

  fit = lmFit(normal_matrix, design)
  fit = contrasts.fit(fit, contrast)
  fit = eBayes(fit)

  result = topTable(fit, adjust.method = "BH", number = Inf)
  result = result |>
    as_tibble(rownames = 'symbol') |>
    rename(padj = adj.P.Val)

  return(result)

}

#' Differential Expression Analysis using DESeq2
#'
#' This function performs differential expression analysis using the DESeq2 package.
#'
#' @param count_matrix A numeric matrix of raw count data, where rows are genes and columns are samples.
#' @param group A factor vector indicating the group assignments for the samples.
#'
#' @return A tibble containing the differential expression results with columns:
#'   - `symbol`: Gene identifiers.
#'   - `logFC`: Log2 fold change.
#'   - `lfcSE`: Standard error of the log2 fold change.
#'   - `stat`: Wald statistic.
#'   - `p.value`: Raw p-value.
#'   - `padj`: Adjusted p-value (FDR).
#'
#' @export
diff_deseq2 = function(count_matrix, group) {

  checkGroup(group)

  library(DESeq2)

  colData = data.frame(group = factor(group))
  dds = DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = colData,
    design = ~ group
  )

  dds = DESeq(dds)
  res = results(dds)

  res_df = as.data.frame(res)
  res_df = res_df[order(res_df$padj), ]

  result = res_df |>
    as_tibble(rownames = 'symbol') |>
    rename(logFC = log2FoldChange)

  return(result)

}

#' Differential Expression Analysis using EdgeR
#'
#' This function performs differential expression analysis using the EdgeR package.
#'
#' @param count_matrix A numeric matrix of raw count data, where rows are genes and columns are samples.
#' @param group A factor vector indicating the group assignments for the samples.
#'
#' @return A tibble containing the differential expression results with columns:
#'   - `symbol`: Gene identifiers.
#'   - `logFC`: Log2 fold change.
#'   - `logCPM`: Log counts per million.
#'   - `LR`: Likelihood ratio statistic.
#'   - `p.value`: Raw p-value.
#'   - `padj`: Adjusted p-value (FDR).
#'
#' @export
diff_edger = function(count_matrix, group) {

  checkGroup(group)

  library(edgeR)

  group = factor(group)
  y = DGEList(counts = count_matrix, group = group)

  keep = filterByExpr(y)
  y = y[keep, , keep.lib.sizes = FALSE]

  y = calcNormFactors(y)
  design = model.matrix(~ group)
  y = estimateDisp(y, design)

  fit = glmFit(y, design)
  lrt = glmLRT(fit, coef = 2)
  res = topTags(lrt, n = Inf)

  res_df = as.data.frame(res)
  result = res_df |>
    as_tibble(rownames = 'symbol') |>
    rename(padj = FDR)

  return(result)

}

