# K-nearest neighbor smoothing based on sNN from clustering
# Inspired by knn-smooth.R: https://github.com/yanailab/knn-smoothing/blob/master/knn_smooth.R
# by Tao Liu <vladimir.liu@gmail.com>

#' kNN-smoothing on a given matrix based on a pre-calculated sNN matrix
#'
#' Note: sparse matrix will be converted to regular matrix, so memory
#' usage will be high.
#'
#' @param M A numeric matrix that will be smoothed. It should have
#'     feature names on rows and cell/barcode names on columns.
#' @param SNN A sNN matrix from graph-based clustering. For example,
#'     the @graph$xxx_sNN matrix from Seurat object. The value can be
#'     regarded as 'similarity' so 1-sNN can be regarded as the
#'     'distance' measurement.
#' @param k Number of nearest neighbours used for smoothing. In
#'     reality, top k+1 nearest neighbors will be considerred
#'     including itself. Default: 10
#' @param method Either 'ave' or 'sum'. The smoothing result will be
#'     either the average of the k-NNs or the sum of the k-NNs. Check
#'     parameter 'strict' as well. Default: 'ave'
#' @param strict A boolean value. If TRUE, those nearest neighbors
#'     with sNN=0 will be excluded, therefore, there may be less than
#'     k NN considerred while calculating the average or sum; if
#'     FALSE, all top nearest neighbors will be included. If 'method'
#'     is 'ave', I recommend 'TRUE' ('FALSE' may be fine too);
#'     however, if 'method' is 'sum', dropping sNN=0 neighbors may
#'     make the cells/barcodes incomparable, so 'FALSE' is recommended
#'     in this case. Default: TRUE
#' @return A smoothed numeric matrix.
#' @export

knn_smoothing <- function(M, SNN, k = 10, method = "ave", strict = TRUE) {
    mat <- as.matrix( M ) # Give up SparseMatrix here.
    cname <- colnames( mat )
    gname <- rownames( mat )
    S <- mat
    if ( method == "ave" ) S <- .smootherAverageNN( mat, SNN, k, strict )
    else S <- .smootherSumNN( mat, SNN, k, strict )
    if (! is.null(cname)) colnames(S) <- cname
    if (! is.null(gname)) rownames(S) <- gname
    return ( S )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param mat A matrix in a shape of features x barcodes.
#' @param SNN A predefined sNN graph in a shape of barcodes x
#'     barcodes.  This SNN matrix can be sparseMatrix, and '.' will be
#'     interpreted as 0. The actual distance between two barcodes is
#'     decided as (1-sNN)
#' @param k An integer to choose k-nearest barcodes (self-exclusive),
#'     so in fact, the top k+1 neighbors including self will be
#'     considered.
#' @param strict A boolean value. If TRUE, those neighbors with sNN=0
#'     will be excluded. Therefore, there may be less than k NN
#'     considerred while calculating the average. Default: TRUE
.smootherAverageNN <- function(mat, SNN, k, strict=TRUE){
    n <- k+1
    sapply(seq_len(ncol(mat)), function(cid){
        nbcid <- head(order(1 - SNN[cid, ]), n) # note: I use 1-sNN as distance
        if ( strict ) nbcid <- nbcid[ SNN[cid, nbcid] >0 ]
        return(Matrix::rowMeans( mat[, nbcid, drop=FALSE] ))
    })
}

#' @param mat A matrix in a shape of features x barcodes.
#' @param SNN A predefined sNN graph in a shape of barcodes x
#'     barcodes.  This SNN matrix can be sparseMatrix, and '.' will be
#'     interpreted as 0. The actual distance between two barcodes is
#'     decided as (1-sNN)
#' @param k An integer to choose k-nearest barcodes (self-exclusive),
#'     so in fact, the top k+1 neighbors including self will be
#'     considered.
#' @param strict A boolean value. If TRUE, those neighbors with sNN=0
#'     will be excluded. Therefore, there may be less than k NN
#'     considerred while computing sum, which is problematic since the
#'     value can not be compared across different barcodes. So the
#'     default is FALSE, and to set it as TRUE is not recommended.
.smootherSumNN <- function(mat, SNN, k, strict=FALSE){
    n <- k+1
    sapply(seq_len(ncol(mat)), function(cid){
        nbcid <- head(order(1 - SNN[cid, ]), n) # note: I use 1-sNN as distance
        if ( strict ) nbcid <- nbcid[ SNN[cid, nbcid] >0 ]
        return(Matrix::rowSums( mat[, nbcid, drop=FALSE] ))
    })
}
