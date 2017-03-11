#' @title Get the OS from R
#' @description Get the type of the operating system. The returned value is used in configuring
#' parallel computation for the implemented algorithms.
#' @export
#' @references  This function is adapted from \url{https://www.r-bloggers.com/identifying-the-os-from-r/}
get_os <- function() {
  # this piece of code is from
  # https://www.r-bloggers.com/identifying-the-os-from-r/
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else {
    ## mystery machine
    # os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
    if (grepl("windows", R.version$os))
      os <- "windows"
  }
  tolower(os)
}



scalar1 <- function(x) {x / sqrt(sum(x^2))}



