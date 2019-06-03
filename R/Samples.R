#' Samples object
#'
#' Samples object elements and methods
#'
#' Creates Samples object, updates its values and contains the functions of
#' this object
#'
#' @param name Name of the object
#' @param data Numeric matrix
#' @param genes name of the rows or genes
#' @param sampleIds name of the columns or samples
#' @param genesToKeep subsample of the genes to keep
#' @param sortGenes .
#' @param name for inner and outerjoin the name of the added Sample object if it was a matrix
#' @param ds another sample object or matrix to do merge to the original sample
#' @import R6
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return a Sample object
#' @examples
#' \dontrun{
#' s <- Samples$new("tcellCD4Profiles", tcellCD4Profiles)
#' s2 <- s$geneSubset(genesToKeep = rownames(tcellCD4Profiles)[1:10])
#' s3 <- s$geneSubset(genesToKeep = rownames(tcellCD4Profiles)[11:20])
#' tpm1 <- s$innerJoin(s2, "s2")
#' colnames(s3$data) <- paste(colnames(s3$data),"s3")
#' tpm2 <- s2$fullOuterJoin(s3)
#' }


Samples <- R6Class("Samples", list(
  # sample name
  name = "sample",
  # columns are samples, rows are genes
  data = NULL,
  # % one column with all gene names.
  #Use name convention "GAPDH" etc and not "ENS..." or "NM..."
  genes = NULL,
  # one row with sample ids
  sampleIds = NULL,
  #usage:
  #l is a logical vector or a vector of indices
  # number of genes
  numberGenes = NULL,
  # number of samples
  numberSamples = NULL,

  initialize = function(name, data = NULL, genes = NULL, sampleIds = NULL, ...) {
    stopifnot(is.character(name), length(name) == 1)
    stopifnot(!is.null(data), is.numeric(data), is.matrix(data))
    self$name <- name
    self$data <- data
    if(is.null(genes)){
      self$genes <- rownames(data)
    } else {
      if(is.character(genes) & length(genes) == dim(data)[1]){
        self$genes <- genes
      } else {
        self$genes <- 1:dim(data)[1]
          warning("The length of the genes array does not match the number of rows in data and/or \n
                it is not character. The genes will be replaced by the name of the rows in data.\n")
      }
    }
    if(is.null(sampleIds)){
      if(is.null(colnames(data))){
        self$sampleIds <- paste(name, 1:dim(data)[2])
      } else{
        self$sampleIds <- colnames(data)
      }
    } else {
      if(is.character(sampleIds) & length(sampleIds) == dim(data)[2]){
        self$sampleIds <- sampleIds
      } else {
        if(is.null(colnames(data))){
          self$sampleIds <- paste(name, 1:dim(data)[2])
        } else{
          self$sampleIds <- colnames(data)
        }
      }
    }
    self$numberSamples <- dim(data)[2]
    self$numberGenes <- dim(data)[1]
    invisible(self)
  },

  print = function(...){
    cat("Sample: \n")
    cat("  name: ", self$name, "\n", sep = "")
    cat("  data:  ", self$data[1,1], "... \n", sep = "")
    cat("  genes:  ", self$genes[1], "... \n", sep = "")
    cat("  number of genes:  ", self$numberGenes, "\n", sep = "")
    cat("  samplesId:  ", self$sampleIds[1:2], "... \n", sep = "")
    cat("  number of samples:  ", self$numberSamples, "\n", sep = "")
    invisible(self)
  },

  #samplesToKeep should be an array, logical array, numeric array or names of samples
  sampleSubset = function(samplesToKeep = NULL, ...){
    stopifnot(is.numeric(samplesToKeep) | is.logical(samplesToKeep) | is.character(samplesToKeep))
    if(length(samplesToKeep) > self$numberSamples){
      warning("The length of the desired samples is larger than the number of samples in the object")
    }
    if(is.character(samplesToKeep)){
      if(sum(!(samplesToKeep %in% self$sampleIds)) > 0){
        warning("Could not find all the samples wanted in the Sample's columns")

      }
    }
    if(is.numeric(samplesToKeep)){
      if(self$numberSamples < max(samplesToKeep) | min(samplesToKeep) <= 0){
        stop("Invalid selection of columns")
      }
    }
    new_sample <- self$sampleIds[self$sampleIds %in% samplesToKeep]
    sNew <- Samples$new(name = paste("subSample", self$name), data =  self$data[,new_sample],
                        genes = rownames(self$data), sampleIds = colnames(self$data[,new_sample]))
    invisible(self)
    return(sNew)
  },

  #genesToKeep should be a vertical cell array, logical array, numeric array or names of genes
  geneSubset = function(genesToKeep = NULL, sortGenes = NULL, ...){
    stopifnot(is.numeric(genesToKeep) | is.logical(genesToKeep) | is.character(genesToKeep))
    if(length(genesToKeep) > self$numberGenes){
      warning("The length of the desired genes is larger than the number of genes in the sample")
    }
    if(is.character(genesToKeep)){
      if(sum(!(genesToKeep %in% self$genes)) > 0){
        warning("Could not find all the genes wanted in the Sample's genes")
      }
    }
    if(is.numeric(genesToKeep)){
      if(self$numberGenes < max(genesToKeep) | min(genesToKeep) <= 0){
        stop("Invalid selection of genes")
      }
    }
    new_genes <- self$genes[self$genes %in% genesToKeep]
    sNew <- Samples$new(name = paste("subSample", self$name), data =  self$data[new_genes,],
                        genes = rownames(self$data[new_genes,]), sampleIds = self$sampleIds)
    invisible(self)
    return(sNew)
  },

  #removes any genes that do not exist in both datasets and scales the columns to sum 1e6
  #returns a Samples object

  innerJoin = function(ds, name = "", scale = TRUE){
    if(!is.matrix(ds) & class(ds)[1] != "Samples"){
      stop("The ds is not a matrix neither a Sample object")
    }
    if(is.matrix(ds)){
      newS <- Samples$new(name, ds)
    } else {
      newS <- ds$clone()
    }
    id <- intersect(self$genes, newS$genes)
    if(length(id) <= 0){
      stop("No overlap between the two datasets")
    } else {
      newData <- cbind(self$data[id,], newS$data[id,])
      if(scale) calcTPM <- tpmDSAVE(newData) else calcTPM <- newData
    }
   return(Samples$new(paste("Intersection", self$name, "and", newS$name), calcTPM))
  },

  #keeps all genes that exist in any dataset and sets them to zero for
  #cells where there is no data
  fullOuterJoin = function(ds, name = "", scale = TRUE){
    if(!is.matrix(ds) & class(ds)[1] != "Samples"){
      stop("The ds is not a matrix neither a Sample object")
    }
    if(is.matrix(ds)){
      newS <- Samples$new(name, ds)
    } else {
      newS <- ds$clone()
    }
    id <- unique(union(rownames(self$data), rownames(newS$data)))
    if(sum(colnames(self$data) %in% colnames(newS$data)) > 0){
      stop("Rename the columns of one of your samples")
    } else {
      idcol <- c(colnames(self$data), colnames(newS$data))
      new_matrix <- matrix(rep(0, length(id) * length(idcol)), ncol = length(idcol))
      rownames(new_matrix) <- id
      colnames(new_matrix) <- idcol
      new_matrix[rownames(self$data), colnames(self$data)] <- self$data[rownames(self$data), colnames(self$data)]
      new_matrix[rownames(newS$data), colnames(newS$data)] <- newS$data[rownames(newS$data), colnames(newS$data)]
      if(scale) calcTPM <- tpmDSAVE(newData) else calcTPM <- newData
      return(Samples$new(paste("Union", self$name, "and", newS$name), calcTPM))
    }
  },

  reScale = function(){
    self$data <- tpmDSAVE(self$data)
  }
))

