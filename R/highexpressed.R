#' highexpressed
#' @aliases high expressed genes libraries 
#' @description Return the highest expressed genes from libraries or treatments.
#'
#' @param obj list object from the function ReadData
#' @param count_zeros boolean of if the number of zero for each library must be shown. Default: TRUE
#' @param num_max int of maximum number of genes reported on high expressed genes. To show all genes expressed use Inf
#' 
#' @return list of the high expressed genes 
#'
#' @keywords High variance high expressed genes  
#'
#' @export
#' 
#' @examples
#' not run
#' dir <-system.file("extdata", package="mnmer")
#' dat <- ReadData (counts, expdata, design="~0+trats")
#' tgs <- highexpressed (dat)
#' 
highexpressed <- function (obj, count_zeros= TRUE, num_max=20)
{
    ctab <- obj$norm_table
    lres <- list()

    gens <- apply(ctab, 1, sum)
    lres$sum_rows <- gens[order(gens, decreasing=TRUE)][1:num_max]

    if (count_zeros)
        lres$number_of_zeros <- apply (ctab, 2, function (x) return (sum(x==0))) 

    tdat <- obj$treat_table

    ltmp <- list()

    for (i in 1:ncol(tdat))
        ltmp[[i]] <- tdat[order(tdat[,i],decreasing = TRUE),i][1:num_max]
    
    names (ltmp) <- colnames(md)
    lres$high_expressed_genes_by_treatment <- ltmp

    ltmp <- list()
    ctab <- obj$norm_table
 
    for (i in 1:ncol(ctab))
        ltmp[[i]] <- ctab[order(ctab[,i], decreasing = TRUE),i][1:num_max]
        
    names (ltmp) <- colnames (ctab)
    lres$high_expressed_genes_by_sample <- ltmp

    return (lres)
}