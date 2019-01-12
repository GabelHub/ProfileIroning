#' Create Results Directories
#'
#' Creates the directories, in which the results are going to be stored.
#' @param homedir The directory in which the result folders will be created.
#' @details All files will be stored in the main folder "Profile-Results". It contains the following subdirectories:
#' \describe{
#'   \item{Figures}{Contains plots of the profiles.}
#'   \item{Fits}{Contains the fitted parameter values of each of the tested models.}
#'   \item{Tables}{Contains the profile likelihood data frames.}
#'   }
#'
#' @export
#' @return Creates storage directories.
#'
#' @examples
#' create.directories(homedir = getwd())

create.directories <- function(homedir) {
  dir <- paste0(homedir, "/Profile-Results")
  dir.create(dir, showWarnings = F)
  dir.create(paste0(dir, "/Figures"), showWarnings = F)
  dir.create(paste0(dir, "/Fits"), showWarnings = F)
  dir.create(paste0(dir, "/Tables"), showWarnings = F)
}
