
# GEO ---------------------------------------------------------------------

#' Read and Process a GEO Series Matrix File
#'
#' This function reads a GEO Series Matrix file from a specified directory, processes the expression
#' and phenotype data, and returns them in tidy formats.
#'
#' @param dir_ A string specifying the directory containing the GEO Series Matrix file. The directory must contain exactly one file.
#' @param alwaysGZ_ A boolean value indicating whether to always decompress gz-compressed txt files in the folder.
#'
#' @return A list with two components:
#'   - `expr`: A tidy tibble of expression data, where rows are samples, and columns are gene symbols.
#'   - `phen`: A tidy tibble of phenotype data, including metadata for each sample.
#'
#' @importFrom GEOquery getGEO
#'
#' @details
#' - The function checks for and processes the platform annotation (GPL). If the annotation cannot be found,
#'   a warning is issued, and the gene IDs are retained as they appear in the expression data.
#' - Genes with the identifier `'permuted_negative'` are removed.
#'
#' @export
readSeries = function(dir_, alwaysGZ_ = T) {

  if (alwaysGZ_) gzDir(dir_)
  series_ = list.files(dir_, full.names = T)

  assert_that(length(series_) == 1)

  gse_ = getGEO(filename = series_, getGPL = F)

  expr_ = gse_ |> Biobase::exprs() |> as_tibble(rownames = 'ID')
  phen_ = pData(gse_) |> as_tibble(rownames = 'sample')

  if (checkGPL(gse_@annotation)) {

    table_id_ = idmap(gse_@annotation, type = 'soft', destdir = 'tmp') |>
      pull(symbol, ID)

  } else {

    warning('GPL Not Found')
    ids_ = unique(expr_$ID)
    table_id_ = setNames(ids_, ids_)

  }

  expr_tidy_ = expr_ |>
    mutate(symbol = table_id_[ID]) |>
    filter(symbol != 'permuted_negative') |>
    distinct(symbol, .keep_all = T) |>
    column_to_rownames('symbol') |>
    select(-ID) |>
    t() |>
    as_tibble(rownames = 'sample')

  phen_tidy_ = expr_tidy_[, 1] |>
    left_join(phen_) |>
    select(sample, title, contains(':ch')) |>
    rename_with(~ .x |> str_sub(end = -5), contains(':ch'))

  return(list(expr = expr_tidy_, phen = phen_tidy_))

}

#' Write GEO Series Data to Disk
#'
#' This function writes processed GEO Series data (expression and phenotype) to specified files in a given directory.
#'
#' @param series_ A list containing two components:
#'   - `expr`: A tibble of expression data.
#'   - `phen`: A tibble of phenotype data.
#' @param dir_ A string specifying the directory where the files will be written. If the directory does not exist, it will be created.
#' @param RData A logical value. If `TRUE`, the function saves the data as an RData file in addition to CSV files. Defaults to `TRUE`.
#' @param gse_ (Optional) A string specifying a GEO Series accession ID (e.g., "GSE12345"). If provided, it is appended to the output filenames.
#'
#' @return None. The function writes files to the specified directory.
#'
#' @details
#' - The function creates the output directory if it does not exist.
#' - Outputs include:
#'   - `expr.csv`: Expression data in CSV format.
#'   - `phen.csv`: Phenotype data in CSV format.
#'   - `all.RData`: An RData file containing the entire dataset (optional).
#'
#' @export
writeSeries = function(series_, dir_, RData = T, gse_ = NULL) {

  if (!dir.exists(dir_)) dir.create(dir_, recursive = T)

  if (!is.null(gse_)) names(series_) = paste(names(series_), gse_, sep = '_')

  cat('Writing expr...\n')
  write_csv(series_$expr, normalizePath(file.path(dir_, 'expr.csv'), '/', F))
  cat('Writing phen...\n')
  write_csv(series_$phen, normalizePath(file.path(dir_, 'phen.csv'), '/', F))

  if (RData) {
    cat('Writing RData...\n')
    env_ = list2env(series_)
    save(list = ls(envir = env_), file = normalizePath(file.path(dir_, 'all.RData'), '/', F), envir = env_)
  }

}
