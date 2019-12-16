#' DSAVEPlotDSAVEScore
#'
#' Plots DSAVE scores, optionally with scores from other datasets as reference.
#'
#' Calculates the average pairwise total variation between two pools of cells with specific size,
#' with a TPM filtration on the genes.
#'
#' @param dsScores Output from DSAVEGetTotalVariation or a list of such outputs
#' @param dsNames The names of the dataset(s), either a string or a list of strings
#' @param ref (Optional) Reference datasets, for example DSAVE::datasetScoresHuman (2000),
#' DSAVE::datasetScoresHuman1000 and DSAVE::datasetScoresHuman500. The reference must have been
#' calculated using the same template as the samples, i.e. match the number of cells.
#' @importFrom methods is as
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return nothing
#' @examples
#' \dontrun{DSAVEPlotDSAVEScore(dsScores, dsNames, DSAVE::datasetScoresHuman)}
#'

DSAVEPlotDSAVEScore <- function(dsScores, dsNames, ref = NA){
  #detect if in param is a list - since the dsData is a list in itself, check for the Rs column...
  if (!is.list(dsScores)) {
    dsScores = list(dsScores)
  }
  if (!is.list(dsNames)) {
    dsNames = list(dsNames)
  }
  if (length(dsScores) != length(dsNames)) {
    stop("Mismatch in length between in params");
  }


  scores = unlist(dsScores)
  names = unlist(dsNames)
  col = rep(0,length(scores))

  if (!is.na(ref[1])) {
    scoresRef <- sapply(ref, function(x) x[[2]])
    namesRef <- sapply(ref, function(x) x[[1]])
    scores = c(scores,scoresRef)
    names = c(names,namesRef)
    col = c(col, rep(1,length(scoresRef)))
  }


  df = data.frame(scores = scores, names=names, col = as.factor(col), stringsAsFactors = F);

  ggplot(df, aes(x = names, fill = col)) +
    geom_bar(aes(weight = scores), show.legend=F) +
    coord_flip() +
    ggtitle("DSAVE BTM Variation Score") +
    ylab("DSAVE BTM score") + xlab("") +
    theme(axis.text.y = element_text(size = 10))
}
