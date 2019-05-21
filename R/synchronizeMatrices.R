#' synchronizeMatrices
#'
#' Not so optimal merge of vectors or matrices
#'
#' Merges vectors or matrices, or vectors and matrices
#'
#' @param mat1 A matrix or vector
#' @param mat2 A matrix or vector
#' @param discardRows a logic vector, TRUE for intersect, FALSE for union
#' @return A matrix or vector
#' @author Juan Inda, <inda@@chalmers.se>
synchronizeMatrices <- function(mat1 = NULL, mat2 = NULL, discardRows = TRUE){
  stopifnot((is.matrix(mat1) & is.numeric(mat1))| is.vector(mat1))
  stopifnot((is.matrix(mat2) & is.numeric(mat2))| is.vector(mat2))

  if(is.vector(mat1) & is.vector(mat2)){
    if(discardRows){
      id <- intersect(mat1, mat2)
      if(length(id) == 0) stop("No intersection between the datasets") else return(id)
    } else {
      return(union(mat1,mat2))
    }
  }
  if(is.vector(mat1) | is.vector(mat2)){
    if(is.vector(mat1) & is.matrix(mat2)){
      M <- mat2
      v <- mat1
    } else {
      if(is.vector(mat2) & is.matrix(mat1)){
        M <- mat1
        v <- mat2
      }
    }
    if(is.null(rownames(M))) rownames(M) <- 1:dim(M)[1]
    if(discardRows){
      id <- intersect(rownames(M), v)
      if(length(id) == 0) stop("No intersection between the datasets") else return(M[id,])
    } else {
      M0 <- matrix(0, ncol = dim(M)[2], nrow = sum(!(v %in% rownames(M))))
      rownames(M0) <- v[!(v %in% rownames(M))]
      return(rbind(M,M0))
    }
  } else {
    if(is.null(rownames(mat1))){
      warning("rownames matrix 1 not setup")
      rownames(mat1) <- 1:dim(mat1)[1]
    }
    if(is.null(rownames(mat2))){
      warning("rownames matrix 1 not setup")
      rownames(mat2) <- 1:dim(mat2)[1]
    }

    if(discardRows){
      id <- intersect(rownames(mat1), rownames(mat1))
      if(length(intersect(colnames(mat1), colnames(mat2)))>0){
        colid <- c(paste(colnames(mat1),"-1"), paste(colnames(mat1),"-2"))
      }  else {
        colid <- c(colnames(mat1), colnames(mat2))
      }
      matRes <- matrix(0, ncol = length(colid), nrow = length(id))
      rownames(matRes) <- id
      colnames(matRes) <- colid
      matRes[id,1:dim(mat1)[2]] <- mat1[id,]
      matRes[id,(dim(mat1)[2] + 1):(dim(mat1)[2] + dim(mat2)[2])] <- mat1[id,]
      return(matRes)
    } else {
      id <- union(rownames(mat1), rownames(mat1))
      if(length(intersect(colnames(mat1), colnames(mat2)))>0){
        colid <- c(paste(colnames(mat1),"-1"), paste(colnames(mat1),"-2"))
      }  else {
        colid <- c(colnames(mat1), colnames(mat2))
      }
      matRes <- matrix(0, ncol = length(colid), nrow = length(id))
      rownames(matRes) <- id
      colnames(matRes) <- colid
      matRes[id %in% rownames(mat1), 1:dim(mat1)[2]] <- mat1[id[id %in% rownames(mat1)],]
      matRes[id %in% rownames(mat2), (dim(mat1)[2] + 1):(dim(mat1)[2]+ dim(mat2)[2])] <- mat2[id[id %in% rownames(mat2)],]
      return(matRes)
    }
  }
}
