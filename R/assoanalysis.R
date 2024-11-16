#' assoanalysis
#' @aliases association analysis  gene expression 
#' @description Perform association rule mining for RNAseq experiments using apriori algorithm
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
assoanalysis <- function (obj, types="apriori")
{
    asso <- list()

    cdat <- obj$count_table


    return (asso)
}