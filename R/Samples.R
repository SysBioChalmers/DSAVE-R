#s <- Samples$new("tcellCD4Profiles", tcellCD4Profiles)
#s2 <- s$geneSubset(genesToKeep = rownames(tcellCD4Profiles)[1:10])
#s3 <- s$geneSubset(genesToKeep = rownames(tcellCD4Profiles)[11:20])
#tpm1 <- s$innerJoin(s2, "s2")

#colnames(s3$data) <- paste(colnames(s3$data),"s3")
#tpm2 <- s2$fullOuterJoin(s3)



Samples <- R6Class("Samples", list(
  # sample name
  name = NULL,
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
  
  initialize = function(name, data = NULL, genes = NULL, sampleIds = NULL) {
    stopifnot(is.character(name), length(name) == 1)
    stopifnot(!is.null(data), is.numeric(data))
    self$name <- name
    self$data <- data
    if(is.null(genes)){
      self$genes <- rownames(data)
    } else {
      if(length(genes) == dim(data)[1]){
        self$genes <- genes
      } else {
        warning("The length of the genes array does not match the number of rows in data.\n
                The genes will be replaced by the name of the rows in data.")
      }
    }
    if(is.null(colnames(data))){
      self$sampleIds <- paste(name, 1:dim(data)[2])
    } else{
      self$sampleIds <- colnames(data)
    }
    
    self$numberGenes <- dim(data)[1]
    invisible(self)
  },
  
  print = function(...){
    cat("Sample: \n")
    cat("  name: ", self$name, "\n", sep = "")
    cat("  data:  ", self$data[1,1:5], "... \n", sep = "")
    cat("  genes:  ", self$genes[1:5], "... \n", sep = "")
    cat("  number of genes:  ", self$numberGenes, "\n", sep = "")
    cat("  samplesId:  ", self$sampleIds[1:2], "... \n", sep = "")
    cat("  number of samples:  ", dim(self$data)[2], "\n", sep = "")
    invisible(self)
  },
  
  #genesToKeep should be a vertical cell array, logical array, numeric array or names of genes
  geneSubset = function(genesToKeep = NULL, sortGenes = NULL, ...){
    stopifnot(is.numeric(genesToKeep) | is.logical(genesToKeep) | is.character(genesToKeep))
    if(length(genesToKeep) > self$numberGenes){
      stop("The length of the desired genes is larger than the number of genes in the sample")
    }
    if(is.character(genesToKeep)){
      if(sum(!(genesToKeep %in% self$genes)) > 0){
        stop("Could not find all the genes wanted in the Sample's genes")
      }
    }
    if(is.numeric(genesToKeep)){
      if(self$numberGenes < max(genesToKeep) | min(genesToKeep) <= 0){
        stop("Invalid selection of genes")
      }
    }
    sNew <- self$clone()
    sNew$data <- sNew$data[genesToKeep,]
    sNew$genes <- rownames(sNew$data)
    sNew$numberGenes <- length(sNew$genes)
    invisible(self)
    return(sNew)
  },
  
  #removes any genes that do not exist in both datasets and scales the columns to sum 1e6
  #returns a Samples object
  
  innerJoin = function(ds, name = ""){
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
      calcTPM <- tpmDSAVE(newData)
    }
   return(Samples$new(paste("Intersection", self$name, "and", newS$name), calcTPM))
  },
  
  #keeps all genes that exist in any dataset and sets them to zero for
  #cells where there is no data
  fullOuterJoin = function(ds, name = ""){
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
      calcTPM <- tpmDSAVE(new_matrix)
      return(Samples$new(paste("Union", self$name, "and", newS$name), calcTPM))
    }
  }
))

