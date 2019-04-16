
#' Download You Too sample data
#' @param folder where downloaded data will be stored
#' @return a list of class \code{brain}
#' @export
#' @examples
#' youtoo <- download_youtoo_data('folder')
#' plot2d(youtoo[[1]], title=name)
#input: folder - location where one wants to store the downloaded data of You Too
#output: list of You Too data of class Brain
download_youtoo_data <- function(folder) {
  youtoo_samples <- c("01", "02", "03", "04", "05", "06", "112", "19", "21", "24")



  dir.create(folder)
  youtoo <- list()
  for(i in (1:length(youtoo_samples))){
    #youtoo
    url <- paste("https://s3.us-east-2.amazonaws.com/deltascope/AT_",youtoo_samples[[i]],"_Probabilities.h5", sep="")
    filePath = paste("./",folder,"/AT_", youtoo_samples[[i]], "_Probabilites.h5", sep="")
    filePath
    download.file(url ,destfile = filePath)
    currFile <- read_h5(filePath)
    youtoo[[i]] <- currFile
  }
  return(youtoo)
}

#' @rdname download_youtoo_data
#' @export

download_wildtype_data <- function(folder) {
  wildtype_samples <- c("101", "102", "103", "104", "105", "106", "107", "108", "109", "110")
  dir.create(folder)
  wildtype <- list()
  for(i in (1:length(wildtype_samples))){
    #wildtype
    url <- paste("https://s3.us-east-2.amazonaws.com/deltascope/AT_wildtype_",wildtype_samples[[i]],"_Probabilities.h5", sep="")
    filePath = paste("./",folder,"/AT_wildtype", wildtype_samples[[i]], "_Probabilites.h5", sep="")
    filePath
    download.file(url ,destfile = filePath)
    currFile <- read_h5(filePath)
    wildtype[[i]] <- currFile
  }
  return(wildtype)
}
