#' ReadData
#' @aliases count matrix read design
#' @description Read data from count matrix and experimental design file.
#'
#' @param counts matrix type containing the gene as rows names and the libraries names as colunms
#' @param col_data data.frame containg the design table of the experiments
#' @param design formula of the factors and treatments
#' @param norm character of the normalization methods: scale (default), TMM, DESeq2
#' @param min_counts minimal row sum that is accpeted by further analisis. Default: 0
#'
#' @return list of objects related to 
#'
#' @keywords High variance replicates formulae
#'
#' @export
#' 
#' @examples
#' not run
#' dir <-system.file("extdata", package="mnmer")
#' dat <- ReadData (counts, expdata, design="~trats")
#' 
ReadData <- function (counts, col_data, design, norm="scale", min_counts=0)
{
	if (!is.matrix(counts)) stop ("The count table must be a matrix type. Please, convert to matrix.")
	if (!is.data.frame(col_data)) stop ("The design must be a data.frame which must contain at leat 2 colunms")
	if (!all(col_data[,1] %in% colnames(counts))) stop ("ERROR between count matrix and the first colunm in the experimental design")

	cdat <- counts[,col_data[,1]]
	rm (counts)

	cdat <- cdat [apply (cdat,1,sum) > min_counts,]
	md <- model.matrix (as.formula(design), data=col_data)

	if (norm == "scale"){
		ff <- apply(cdat,2,sum)
		ff <- ff/min(ff)
		cdat <- t(t(cdat)/ff)
		tdat <- ctab %*% md
    	ff <- apply (md,2,sum)
    	tdat <- t(t(tdat)/ff)
	} else {
       stop ("Normalization other than scale is not implemented yet!")
    }
	

	return (list(norm_table=cdat, treat_table=tdat, model_matrix=md, norm=norm, min_counts=0))
}

