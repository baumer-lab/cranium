#' Download YouToo sample data
#' @param path where downloaded data will be stored
#' @param type \code{wildtype} or \code{youtoo}. Default is \code{wildtype}
#' @param pattern regular expression to match file names
#' @param ... arguments passed to \code{\link{download.file}}
#' @return a list of class \code{brain}
#' @export
#' @examples
#' wildtype_files <- download_brains(pattern = "101")
#' youtoo_files <- download_brains(type = "youtoo", pattern = "6")
#' brains <- read_brains(c(wildtype_files, youtoo_files))
#' plot2d(brains[[1]], title = name)

download_brains <- function(path = ".", type = "wildtype", pattern = NULL, ...) {
  youtoo_samples <- c("01", "02", "03", "04", "05", "06", "112", "19", "21", "24")
  #youtoo
  wildtype_samples <- c("101", "102", "103", "104", "105", "106", "107", "108", "109", "110")

  x <- tibble::tibble(
    src = paste0(
    "https://s3.us-east-2.amazonaws.com/deltascope/AT_",
    if (type == "wildtype") {
      "wildtype_"
    } else {
      ""
    },
    if (type == "wildtype") {
      wildtype_samples
    } else {
      youtoo_samples
    },
    "_Probabilities.h5"),
    lcl = fs::path(path, basename(src)),
    exists = fs::file_exists(lcl)
  )
  if (!is.null(pattern)) {
    x <- x %>%
      filter(grepl(pattern, src))
  }
  fs::dir_create(path)

  message(paste("Downloading", sum(!x$exists), "new files. ",
                sum(x$exists), "untouched."))
  x %>%
    filter(!exists) %>%
    group_split(src) %>%
    purrr::walk(~download.file(.$src, destfile = .$lcl))

  return(fs::dir_ls(path, regexp = pattern))
}

#' @rdname download_brains
#' @param file path to HDF5 file
#' @export

read_brains <- function(file) {
  purrr::map(file, read_h5)
}

