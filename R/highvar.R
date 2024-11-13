#' highvar
#' @aliases high variance cluster PCA gene expression 
#' @description Analyse if a experiments has a high variance among replicates and treatments.
#'
#' @param obj list object from the function ReadData
#' @param types character of type of nalysis to be performed 
#' 
#' @return list of the treatments and its varince to compared
#' @keywords high variance PCA cluster association
#'
#' @export
#' 
#' @examples
#' not run
#' dir <-system.file("extdata", package="mnmer")
#' dat <- ReadData (counts, expdata, design="~trats")
#' tgs <- varanalysis (dat)
#' 
highvar <- function (obj, type="PCA"){
    lhvar <- list()




    return (lhvar)
}