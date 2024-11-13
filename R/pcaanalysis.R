#' pcaanalysis
#' @aliases principal component analysis gene expression 
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
pcaanalysis <- function (obj, plotcomps=TRUE, plotpcacor=TRUE)
{
    lres <- list()

    ctab <- obj$norm_table
    res.pca <- PCA(ctab, graph = FALSE)

    tdat <- obj$treat_table
    res2.pca <- PCA(tdat, graph = FALSE)

    lres$eigenvalue_percentage_of_variance_cumulative_by_samples  <- res.pca$eig

    if (plotcomps){
        p1 <- fviz_eig(res.pca) + ggtitle ("Components by Samples")
        p2 <- fviz_eig(res2.pca) + ggtitle ("Components by Treatments")
        grid.arrange(p1, p2, nrow = 1)
    }

    if (plotpcacor) {
        p3 <- fviz_pca_var(res.pca, col.var = "black") + ggtitle ("Samples")
        p4 <- fviz_pca_var(res2.pca, col.var = "black") + ggtitle ("Treatments")
        grid.arrange(p3, p4, nrow = 1)
    }

    return (lres)

}