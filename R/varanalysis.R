#' varanalysis
#' @aliases variance standard deviation gene expression 
#' @description Return the variance analysis including near zero variance with treatments.
#'
#' @param obj list object from the function ReadData
#' @param num int of number of standard deviation reported. Default: 20. To show all genes expressed use Inf
#' 
#' @return list of the variane by genes among treatments 
#'
#' @keywords High variance high expressed genes  
#'
#' @export
#' 
#' @examples
#' not run
#' dir <-system.file("extdata", package="mnmer")
#' dat <- ReadData (counts, expdata, design="~trats")
#' tgs <- varanalysis (dat)
#' 
varanalysis <- function (obj, num=20){
    ctab <- obj$norm_table
    md <- obj$model_matrix

    lres <- list()

    tmp <- apply (ctab, 1 , sd)
    lres$sd_all_rows <- tmp[order(tmp,decreasing = TRUE)][1:num]

    lres$near_zero_variance <- nearZeroVar (ctab, saveMetrics=TRUE)

    return (lres)
}