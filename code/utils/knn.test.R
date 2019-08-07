knn.test <- function(sample, time, knn, correlate = TRUE) {
  # source("pair.consecutive.R")
  consec <- pair.consecutive(sample, time)
  consec <- mutate(consec,
                   sample.x = as.character(sample.x),
                   sample.y = as.character(sample.y))
  consec <- mutate(consec, delta.t = time.y - time.x)
  names(knn) <- sample
  consec <- mutate(consec,
                   knn.x = knn[sample.x],
                   knn.y = knn[sample.y],
                   delta.knn = knn.y - knn.x)
  if (correlate) {
    cor(consec$knn.x, consec$delta.knn)
  } else {
    consec
  }
}