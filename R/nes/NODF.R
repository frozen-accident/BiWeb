NODF <- function(matrix, nodf_strict = 1) {
  # TO DO, implement the decreasing strict (nodf_strict) section

  m <- 1.0 * matrix
  n_rows <- nrow(m)
  n_cols <- ncol(m)

  C1 <- m %*% t(m)
  D1 <- matrix(diag(C1), n_rows, n_rows)
  nr <- C1 * (D1 != t(D1))
  nr <- nr / pmin(D1, t(D1))

  C2 <- t(m) %*% m
  D2 <- matrix(diag(C2), n_cols, n_cols)
  nc <- C2 * (D2 != t(D2))
  nc <- nc / pmin(D2, t(D2))

  nr <- 100 * sum(nr, na.rm = TRUE) / 2
  nc <- 100 * sum(nc, na.rm = TRUE) / 2
  nvalue <- nr + nc

  denom <- n_rows * (n_rows - 1) / 2 + n_cols * (n_cols - 1) / 2
  nr <- nr / (n_rows * (n_rows - 1) / 2)
  nc <- nc / (n_cols * (n_cols - 1) / 2)
  nvalue <- nvalue / denom

  list(nodf_rows = nr, nodf_cols = nc, nodf = nvalue)
}
