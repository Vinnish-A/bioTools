
# utils -------------------------------------------------------------------


## colors ----

#' @export
color_morandi = c('#8ac3c6', '#999fbf', '#afafad', '#87b8de', '#f38684', '#a48999')
#' @export
color_macaron = c("#E45D61", "#4A9D47", "#F19294", "#96C3D8", "#5F9BBE", '#394c81', "#F5B375", "#67A59B", '#94697a', "#A5D38F", '#999fbf')
#' @export
color_monika = c('#999fbf', '#94697a', '#4A9D47')
#' @export
color_myth = c('#026401', 'black', '#cc0303')

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

#' @export
withCal = with

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
