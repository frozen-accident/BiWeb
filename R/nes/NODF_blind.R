NODF_blind <- function(matrix) {

  matrix <- matrix != 0
  nrows <- nrow(matrix)
  ncols <- ncol(matrix)

  # sort_matrix
  rowsort <- order(-rowSums(matrix))
  colsort <- order(-colSums(matrix))
  matrix <- matrix[rowsort, colsort, drop = FALSE]

  sumrows <- rowSums(matrix)
  sumcols <- colSums(matrix)

  nr <- 0
  # Fill for columns
  for (i in seq_len(nrows)) {
    for (j in colon(i + 1, nrows)) {
      if (sumrows[j] < sumrows[i] && sumrows[j] > 0) {
        nr <- nr + 100 * sum(matrix[i, matrix[j, ]]) / sumrows[j]
      }
    }
  }

  nc <- 0
  for (k in seq_len(ncols)) {
    for (l in colon(k + 1, ncols)) {
      if (sumcols[l] < sumcols[k] && sumcols[l] > 0) {
        nc <- nc + 100 * sum(matrix[matrix[, l], k]) / sumcols[l]
      }
    }
  }

  nvalue <- nr + nc
  denom <- nrows * (nrows - 1) / 2 + ncols * (ncols - 1) / 2
  nr <- nr / (nrows * (nrows - 1) / 2)
  nc <- nc / (ncols * (ncols - 1) / 2)
  nvalue <- nvalue / denom

  imagesc(matrix)

  list(nvalue = nvalue, nr = nr, nc = nc)
}
