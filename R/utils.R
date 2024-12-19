
# utils -------------------------------------------------------------------

myApply = function(X, MARGIN, FUN, ..., simplify = TRUE) {

  if ('future.apply' %in% installed.packages()) {

    if (future::nbrOfWorkers() > 1) {
      res_ = future.apply::future_apply(X, MARGIN, FUN, ..., simplify = TRUE)
    } else {
      res_ = apply(X, MARGIN, FUN, ..., simplify = TRUE)
    }

  } else {

    res_ = apply(X, MARGIN, FUN, ..., simplify = TRUE)

  }

  return(res_)

}

## colors ----

#' Convert Character Vector to Colored Character Class
#'
#' This function converts a character vector to a special class `coloredChr`, which is used to store and print color data.
#'
#' @param vec_ A character vector to be converted.
#' @return A character vector with the class attribute `coloredChr`.
#'
#' @examples
#' vec = c("red", "green", "blue")
#' colored_vec = asColoredChr(vec)
#' class(colored_vec)  # "coloredChr"
#'
#' @export
asColoredChr = function(vec_) {

  if (!is.character(vec_)) stop("Input must be a character vector")
  class(vec_) = c("coloredChr", class(vec_))

  return(vec_)

}

#' Print Colored Character Vector with ANSI Colors
#'
#' This function prints a character vector with ANSI colors in a plot and also displays the vector values in the console.
#' It is designed to work with objects of class `coloredChr`.
#'
#' @param x An object of class `coloredChr`, which is a character vector containing color information.
#' @param ... Additional arguments passed to the method (currently not used).
#' @return No return value. The function prints a plot and the color vector in the console.
#'
#' @method print coloredChr
#'
#' @import ggplot2 dplyr
#'
#' @examples
#' vec = asColoredChr(c("red", "green", "blue"))
#' print(vec)  # This will display a plot and print the vector to the console.
#'
#' @export
print.coloredChr = function(x, ...) {

  cat("Color Vector with ANSI Colors:\n")

  dp_ = data.frame(
    x = 1:length(x),
    y = rep(1, length(x)),
    color = x
  )

  p_ = dp_ |>
    ggplot(aes(x = x, y = y, fill = color)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_identity() +
    theme_void() +
    theme(legend.position = "none")

  print(p_)
  cat(x, sep = "\n")

}

dict_color_1 = lst(
  conceptualization = '#999fbf',
  presentation = '#94697a'
)

dict_color_2 = lst(
  rb = c('red', 'blue'),
  rb_classic = c('firebrick', 'navy'),
  rb_soft = c("#f19294", "#5F9BBE"),
  ddlc = c('#94697a', '#999fbf')
) |> lapply(rev)

dict_color_3 = lst(
  htmap = c('#026401', 'black', '#cc0303'),
  monika = c('#999fbf', '#94697a', '#57a5a4')
)

dict_color_n = lst(
  morandi = c('#8ac3c6', '#999fbf', '#afafad', '#e1d5d7', '#87b8de', '#f38684', '#a48999'),
  macaron = c("#F19294", "#A5D38F", "#96C3D8", "#5F9BBE", "#E45D61", "#4A9D47", '#394c81', "#F5B375", "#67A59B", '#999fbf', '#94697a', '#6a3d9a', '#cab2d6', '#fdbf6f', '#ff7f00', '#ffff99', '#b15928')
)

#' Get Color from dict_color_1, dict_color_2, dict_color_3, and dict_color_n
#'
#' This function retrieves the color from the dictionary based on the specified `which_` parameter.
#' The default is 'conceptualization'.
#'
#' @param which_ A character string specifying the key to retrieve from the color dictionary. Defaults to 'conceptualization'.
#' @return A color corresponding to the specified key in the dictionary.
#'
#' @keywords internal
get_color_1 = function(which_ = 'conceptualization') {

  dict_ = c(dict_color_1, dict_color_2, dict_color_3, dict_color_n)

  dict_[[which_]]

}

#' Get Color from dict_color_2, dict_color_3, and dict_color_n
#'
#' This function retrieves the color from the dictionary based on the specified `which_` parameter.
#' The default is 'ddlc'.
#'
#' @param which_ A character string specifying the key to retrieve from the color dictionary. Defaults to 'ddlc'.
#' @return A color corresponding to the specified key in the dictionary.
#'
#' @keywords internal
get_color_2 = function(which_ = 'ddlc') {

  dict_ = c(dict_color_2, dict_color_3, dict_color_n)

  dict_[[which_]]

}

#' Get Color from dict_color_3 and dict_color_n
#'
#' This function retrieves the color from the dictionary based on the specified `which_` parameter.
#' The default is 'htmap'.
#'
#' @param which_ A character string specifying the key to retrieve from the color dictionary. Defaults to 'htmap'.
#' @return A color corresponding to the specified key in the dictionary.
#'
#' @keywords internal
get_color_3 = function(which_ = 'htmap') {

  dict_ = c(dict_color_3, dict_color_n)

  dict_[[which_]]

}

#' Get Color from dict_color_n
#'
#' This function retrieves the color from the dictionary based on the specified `which_` parameter.
#' The default is 'macaron'.
#'
#' @param which_ A character string specifying the key to retrieve from the color dictionary. Defaults to 'macaron'.
#' @return A color corresponding to the specified key in the dictionary.
#'
#' @keywords internal
get_color_n = function(which_ = 'macaron') {

  dict_ = dict_color_n

  dict_[[which_]]

}

#' Get a specified number of colors from a predefined dictionary of colors
#'
#' This function retrieves `n_` colors from the color dictionary. If `n_` is greater than 3,
#' it retrieves from the `dict_color_n` dictionary; otherwise, it selects from the respective predefined color sets.
#'
#' @param n_ A numeric value specifying the number of colors to retrieve. Must be between 1 and 17.
#' @param which_ A character string specifying which color to retrieve. Defaults to NULL, in which case all colors are retrieved.
#' @return A colored vector of the requested number of colors.
#'
#' @examples
#' get_color(3, 'morandi')
#' get_color(5, 'morandi')
#'
#' @export
get_color = function(n_, which_ = NULL) {

  if (n_ < 1 | n_ > 17) stop('nColor < 1 or nColor > 17')

  if (n_ > 3) {
    suffix_ = 'n'
  } else {
    suffix_ = n_
  }

  if (is.null(which_)) {
    res_ = do.call(paste0('get_color_', suffix_), list())
  } else {
    res_ = do.call(paste0('get_color_', suffix_), list(which_ = which_))
  }

  res_ = asColoredChr(res_[1:n_])

  return(res_)

}

## fun ----

#' Append Elements to a List by Name
#'
#' This function appends new elements to an existing list, using the provided names for the new elements.
#'
#' @param lst_ A list to which new elements will be appended.
#' @param ... Named arguments representing the elements to be appended to the list.
#'
#' @return A list with the new elements appended.
#'
#' @examples
#' my_list = list(a = 1, b = 2)
#' my_list = appendWithName(my_list, c = 3, d = 4)
#' print(my_list)
#'
#' @export
appendWithName = function(lst_, ...) {

  lst_appending_ = list(...)
  for (i_ in seq_along(lst_appending_)) {

    name_ = names(lst_appending_)[[i_]]
    value_ = lst_appending_[[i_]]

    lst_[[name_]] = value_

  }

  return(lst_)

}

#' Execute a Function in the Context of an Object
#'
#' This function evaluates a specified function with additional parameters evaluated in the context of an object.
#'
#' @param x_ An object that provides the context for evaluating the parameters.
#' @param fun_ A function to be executed.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return The result of the function execution.
#'
#' @examples
#' my_env = new.env()
#' my_env$a = 2
#' my_env$b = 3
#' my_env |> doCal(sum, a, b)
#'
#' @export
doCal = function(x_, fun_, ...) {

  params_ = lapply(as.list(match.call())[-c(1:3)], \(x__) eval(x__, envir = x_))
  # if (is.null(names(params_)) | '' %in% names(params_)) names(params_) = names(formals(fun_))
  res_ = do.call(fun_, params_)

  return(res_)

}

#' Evaluate an Expression within a Specific Environment
#'
#' This function evaluates an expression in the environment of a specified data object, which allows for customized scoping and evaluation.
#'
#' @param x_ An environment or data object in which the expression should be evaluated.
#' @param expr_ An R expression to be evaluated within the context of `x_`.
#' @return The result of evaluating the expression `expr_` within the environment `x_`.
#'
#' @examples
#' data.frame(a = 1:3, b = 4:6) |> withCal(a + b)
#'
#' @export
withCal = function(x_, expr_) {

  eval(substitute(expr_), x_, enclos = parent.frame())

}

#' Apply a Function to Specific Elements of a List
#'
#' This function applies a specified function to either a subset of elements (lucky ones) or the rest (unlucky ones) in a list.
#'
#' @param lst_ A list whose elements the function will be applied to.
#' @param fun_ A function to apply to the elements.
#' @param luckyOnes_ An integer vector indicating the indices of the lucky elements to apply the function to.
#' @param reverse_ A logical value. If `TRUE`, the function is applied to the unlucky elements instead. Defaults to `FALSE`.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A list with the function applied to the specified elements.
#'
#' @examples
#' my_list = list(1, 2, 3, 4)
#' biasedMap(my_list, \(x) x^2, luckyOnes_ = c(2, 3))
#'
#' @export
biasedMap = function(lst_, fun_, luckyOnes_ = 1, reverse_ = F, ...) {

  luckyOnes_[luckyOnes_ < 0] = luckyOnes_[luckyOnes_ < 0] + 1 + length(lst_)
  luckyOnes_ = sort(luckyOnes_)
  unluckyOnes_ = setdiff(1:length(lst_), luckyOnes_)

  if (!reverse_) {
    lst_[luckyOnes_] = lapply(lst_[luckyOnes_], fun_, ...)
  } else {
    lst_[unluckyOnes_] = lapply(lst_[unluckyOnes_], fun_, ...)
  }

  return(lst_)

}

#' Extract and Uncompress Gzipped Files in a Directory
#'
#' This function uncompresses all `.gz` files in a specified directory.
#'
#' @param dir_ A string specifying the directory containing `.gz` files.
#'
#' @return Invisible `NULL`.
#'
#' @details
#' - Removes the original `.gz` files after extraction.
#' - If no `.gz` files are found, prints "NULL" and does nothing.
#'
#' @examples
#' gzDir("path/to/directory")
#'
#' @export
gzDir = function(dir_) {

  filenames_ = list.files(dir_, full.names = T, pattern = 'gz')
  if (length(filenames_) == 0) {
    cat('NULL')
    invisible(NULL)
  } else {
    lapply(filenames_, R.utils::gunzip, remove = T)
    invisible(NULL)
  }

}

#' Safely Load Data from an RData File
#'
#' This function loads objects from an RData file into a list, avoiding pollution of the global environment.
#'
#' @param path2RData_ A string specifying the path to the RData file.
#'
#' @return A named list containing all objects from the RData file.
#'
#' @examples
#' # Save and load an example RData file
#' my_data = list(a = 1, b = 2)
#' save(my_data, file = "example.RData")
#' loaded_data = loadSafely("example.RData")
#' print(loaded_data)
#'
#' @export
loadSafely = function(path2RData_) {

  env_ = environment()
  load(path2RData_, envir = env_)
  name_ = setdiff(ls(envir = env_), c('env_', 'path2RData_'))
  lst_ = mget(name_, envir = env_)

  return(lst_)

}


