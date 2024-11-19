
# tidy --------------------------------------------------------------------


#' Calculate Row Sums and Add as a New Column
#'
#' This function calculates the row sums of selected numeric columns in a data frame
#' and adds the result as a new column with a specified name.
#'
#' @param df_ A data frame containing the data.
#' @param name_ A string specifying the name of the new column where the row sums will be stored.
#' @param con_ An expr about how to select columns.
#'             Defaults to `where(is.numeric)`, which selects all numeric columns.
#'
#' @return A modified data frame with the new column containing row sums.
#'
#' @import dplyr
#' @import rlang
#'
#' @examples
#' library(dplyr)
#' df = tibble(a = 1:3, b = 4:6, c = letters[1:3])
#' df |> tidyRowSums("row_sum")
#'
#' @export
tidyRowSums = function(df_, name_, con_ = where(is.numeric)) {

  con_ = enquo(con_)

  df_[[name_]] = rowSums(df_ |> select(!!con_))

  return(df_)

}

#' Transpose a tibble
#'
#' This function transposes a data frame, converting its columns into rows,
#' and optionally renames the rownames column in the output tibble.
#'
#' @param df_ A data frame to be transposed.
#' @param id_new_ A string specifying the name of the new column that will contain the original column names (now row names).
#' @param id_old_ A string specifying the name of the column to be used as rownames before transposing.
#'                Defaults to the first column of the input data frame.
#'
#' @return A tibble with the transposed data, including a column with the original column names.
#'
#' @import tibble
#'
#' @examples
#' df = tibble(id = c("row1", "row2"), a = c(1, 2), b = c(3, 4))
#' df |> tidyT("new_id")
#'
#' @export
tidyT = function(df_, id_new_, id_old_ = colnames(df_)[[1]]) {

  df_ |>
    column_to_rownames(id_old_) |>
    t() |>
    as_tibble(rownames = id_new_)

}
