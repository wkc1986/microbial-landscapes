#' Make pairs of consecutive sample IDs and times
#'
#' @param sample
#' @param time
#'
#' @return
#' @export
#'
#' @examples
pair.consecutive <- function(sample, time) {
  df <- data.frame(sample.x = sample, time.x = time)
  df <- df[order(df$time.x),]
  df1 <- df[-nrow(df),]
  df2 <- df[-1,]
  names(df2) <- c("sample.y", "time.y")
  cbind(df1, df2)
}