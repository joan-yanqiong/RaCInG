#' Get name of file without extension
#'
#' @param filepath Path to file
#' @return Name of file without extension
#'
#' @examples get_name("R/my_script.R") returns 'my_script'
#' @importFrom tools file_path_sans_ext
#' @importFrom fs path_file
#' @export
get_name <- function(filepath) {
  return(file_path_sans_ext(path_file(filepath)))
}
#' Create directory if it does not exist
#'
#' @param dir_path Path to directory to be created
#'
#' @examples create_dir("/Users/johndoe/my_new_dir")
#' @export
create_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    log_info(glue("Creating directory {dir_path}"))
    dir.create(dir_path, recursive = TRUE)
  } else {
    log_warn(glue("Directory {dir_path} already exists."))
  }
}

#' Scale values to [0, 1]
#'
#' @param x numeric vector
#' @return vector with values scaled to [0, 1]
#' @export
#'
#' @examples range01(c(1, 2, 3))
range01 <- function(x) {
  return(x - min(x)) / (max(x) - min(x))
}
