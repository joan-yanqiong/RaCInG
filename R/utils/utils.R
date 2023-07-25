# \' Get name of file without extension
# \'
# \' @param filepath
# \' @return Name of file without extension
# \'
# \' @examples get_name("R/my_script.R") returns 'my_script'
get_name <- function(filepath) {
  return(tools::file_path_sans_ext(fs::path_file(filepath)))
}
# \' Create directory if it does not exist
# \'
# \' @param dir_path Path to directory to be created
# \'
# \' @examples create_dir("/Users/johndoe/my_new_dir")
create_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    log_info(glue("Creating directory {dir_path}"))
    dir.create(dir_path, recursive = TRUE)
  } else {
    log_warn(glue("Directory {dir_path} already exists."))
  }
}
# \' Set up logging
# \'
# \' @param log_level level of logging (1-5)
# \' @return log4r logger object
# \'
# \' @examples logr <- init_logging(3)
init_logging <- function(log_level) {
  log_level_options <- c(
    `1` = "FATAL", `2` = "ERROR", `3` = "WARN", `4` = "INFO",
    `5` = "DEBUG"
  )
  return(log4r::logger(
    threshold = log_level_options[as.character(log_level)],
    appenders = log4r::console_appender(layout = log4r::default_log_layout())
  ))
}
# \' Logging functions: log_info
# \'
# \' @param logr
# \'
# \' @examples log_info("Hello world!")
log_info <- function(...) {
  log4r::info(logr, paste0(...))
}

# \' Logging functions: log_error
# \'
# \' @param logr
# \'
# \' @examples log_error("Hello world!")
log_error <- function(...) {
  log4r::error(logr, paste0(...))
}

# \' Logging functions: log_fatal
# \'
# \' @param logr
# \'
# \' @examples log_fatal("Hello world!")
log_fatal <- function(...) {
  log4r::fatal(logr, paste0(...))
}

# \' Logging functions: log_debug
# \'
# \' @param logr
# \'
# \' @examples log_debug("Hello world!")
log_debug <- function(...) {
  log4r::debug(logr, paste0(...))
}

# \' Logging functions: log_warn
# \'
# \' @param logr
# \'
# \' @examples log_warn("Hello world!")
log_warn <- function(...) {
  log4r::warn(logr, paste0(...))
}


pairs_oi <- list(
  c("Macrophage", "AC"), c("Macrophage", "MES"), c("Macrophage", "NPC & OPC")
)
pairs_oi_rev <- lapply(pairs_oi, function(p) {
  return(rev(p))
})
all_pairs_oi <- sapply(c(pairs_oi, pairs_oi_rev), function(p) {
  return(paste(p, collapse = "__"))
})
