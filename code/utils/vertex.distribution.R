#' Group by variables and get vertex distribution
#'
#' @param vertex     vector of vertex IDs
#' @param group.vars DF of grouping variables
#'
#' @return DF of counts + frequencies per vertex + group
#' @export
#'
#' @examples
vertex.distribution <- function(vertex, group.vars) {
 df <- data.frame(vertex = vertex, group.vars)
 grouped <- group_by_all(df)
 smry <- summarise(grouped, n = n())
 smry <- ungroup(smry)
 smry <- group_by_at(smry, names(group.vars))
 smry <- mutate(smry, f = n / sum(n))
 ungroup(smry)
}