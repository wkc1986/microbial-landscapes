#' Make square distance matrix from vectors of sample ids and distance
#'
#' @param sample.x
#' @param sample.y
#' @param d
#'
#' @return
#' @export
#'
#' @examples
dlist2dm <- function(sample.x, sample.y, d) {
  df <- data.frame(sample.x, sample.y, d)
  dm <- dcast(df, sample.x ~ sample.y, value.var = "d")
  rn <- dm$sample.x
  dm <- dm[, -1]
  dm <- data.matrix(dm)
  rownames(dm) <- rn
  dm
}