#' JS divergence between vertex distributions of states
#'
#' @param f      vector of vertex frequencies
#' @param vertex vector of vertex IDs
#' @param state  vector of state labels
#'
#' @return
#' @export
#'
#' @examples
state.divergence <- function(f, vertex, state) {
  library(philentropy)
  vd <- data.frame(f, vertex, state)
  m <- dcast(vd, state ~ vertex, value.var = "f", fill = 0)
  rownames(m) <- m$state
  m <- m[, -1] %>%
    data.matrix
  sd <- JSD(m)
  rownames(sd) <- rownames(m)
  colnames(sd) <- rownames(m)
  melt(sd, varnames = c("state.x", "state.y"), value.name = "JSD")
}