#' Set up logging
#'
#' @param log_level level of logging (1-5)
#' @param log_file path to log file (optional)
#' @return logger object
#'
#' @examples logr <- init_logging(3)
#' @export
#' @importFrom log4r logger console_appender file_appender default_log_layout
init_logging <- function(log_level, log_file = NULL) {
    log_level_options <- c(
        `1` = "FATAL", `2` = "ERROR", `3` = "WARN", `4` = "INFO",
        `5` = "DEBUG"
    )
    if (!is.null(log_file)) {
        console_app <- console_appender(layout = default_log_layout())
        file_app <- file_appender(log_file,
            append = FALSE,
            layout = default_log_layout()
        )
        return(logger(
            threshold = log_level_options[as.character(log_level)],
            appenders = list(console_app, file_app)
        ))
    }
    return(logger(
        threshold = log_level_options[as.character(log_level)],
        appenders = console_appender(layout = default_log_layout())
    ))
}

#' Logging functions: log_info
#'
#' @param logr logger object
#' @param ... message
#'
#' @examples log_info("Hello world!")
#' @export
#' @importFrom log4r info
log_info <- function(...) {
    info(logr, paste0(...))
}

#' Logging functions: log_error
#'
#' @param logr logger object
#' @param ... message
#'
#' @examples log_error("Hello world!")
#'
#' @export
#' @importFrom log4r error
log_error <- function(...) {
    error(logr, paste0(...))
}

#' Logging functions: log_fatal
#'
#' @param logr logger object
#' @param ... message
#'
#' @examples log_fatal("Hello world!")
#'
#' @export
#' @importFrom log4r fatal
log_fatal <- function(...) {
    fatal(logr, paste0(...))
}

#' Logging functions: log_debug
#'
#' @param logr logger object
#' @param ... message

#' @examples log_debug("Hello world!")
#'
#' @export
#' @importFrom log4r debug
log_debug <- function(...) {
    debug(logr, paste0(...))
}

#' Logging functions: log.warn
#'
#' @param logr logging object
#' @param ... message
#'
#' @examples log_warn("Hello world!")
#' @export
#' @importFrom log4r warn
log_warn <- function(...) {
    log4r::warn(logr, paste0(...))
}
