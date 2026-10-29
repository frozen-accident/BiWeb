NullModels <- local({

  NULL_1 <- function(MatrixOrSizeRows, n_cols, ones) {

    if (missing(n_cols)) {
      n_rows <- nrow(MatrixOrSizeRows)
      n_cols <- ncol(MatrixOrSizeRows)
      ones <- sum(MatrixOrSizeRows > 0)
    } else {
      n_rows <- MatrixOrSizeRows
    }

    x <- matrix(0, n_rows, n_cols)
    x[sample.int(n_rows * n_cols, ones)] <- 1
    x
  }

  NULL_2 <- function(n_rows, n_cols, probRows, probCols) {

    if (missing(n_cols)) {
      matrix <- n_rows
      n_rows <- nrow(matrix)
      n_cols <- ncol(matrix)
      probCols <- colSums(matrix) / n_rows
      probRows <- rowSums(matrix) / n_cols
    }

    if (n_rows != length(probRows)) {
      stop("The row size is different from the assigned probabilities")
    }
    if (n_cols != length(probCols)) {
      stop("The column size is different from the assigned probabilities")
    }

    pr <- matrix(probRows, n_rows, n_cols)
    pc <- matrix(probCols, n_rows, n_cols, byrow = TRUE)

    mat <- (pr + pc) / 2

    matrix(runif(n_rows * n_cols), n_rows, n_cols) <= mat
  }

  NULL_3ROW <- function(n_rows, n_cols, probRows) {
    # Return a random network in which rows fill according to probRows

    if (missing(n_cols)) {
      matrix <- n_rows
      n_rows <- nrow(matrix)
      n_cols <- ncol(matrix)
      probRows <- rowSums(matrix) / n_cols
    }

    if (n_rows != length(probRows)) {
      stop("The row size is different from the assigned probabilities")
    }

    repeat {
      matrix <- matrix(runif(n_rows * n_cols), n_rows, n_cols) < matrix(probRows, n_rows, n_cols)
      if (all(colSums(matrix) > 0) && all(rowSums(matrix) > 0)) break
    }
    matrix
  }

  NULL_3COL <- function(n_rows, n_cols, probCols) {

    if (missing(n_cols)) {
      matrix <- n_rows
      n_rows <- nrow(matrix)
      n_cols <- ncol(matrix)
      probCols <- colSums(matrix) / n_rows
    }

    if (n_cols != length(probCols)) {
      stop("The column size is different from the assigned probabilities")
    }

    repeat {
      matrix <- matrix(runif(n_rows * n_cols), n_rows, n_cols) < matrix(probCols, n_rows, n_cols, byrow = TRUE)
      if (all(colSums(matrix) > 0) && all(rowSums(matrix) > 0)) break
    }
    matrix
  }

  # Name of a null model, given either the function itself or its name
  # (e.g. NullModels$NULL_2, "NULL_2" or "NullModels.NULL_2").
  MODEL_NAME <- function(model) {
    if (is.character(model)) return(sub("^NullModels[.$]", "", model))
    for (nm in c("NULL_1", "NULL_2", "NULL_3ROW", "NULL_3COL")) {
      if (identical(model, get(nm))) return(nm)
    }
    "NULL_1"
  }

  NULL_MODEL <- function(adjmatrix, model, replicates) {

    # R's RNG is already randomly seeded per session; call set.seed() for reproducibility.

    adjmatrix <- adjmatrix != 0

    n_edges <- sum(adjmatrix)
    n_rows <- nrow(adjmatrix)
    n_cols <- ncol(adjmatrix)
    p_cols <- colSums(adjmatrix) / n_rows
    p_rows <- rowSums(adjmatrix) / n_cols

    nullmodel <- MODEL_NAME(model)

    rmatrices <- vector("list", replicates)
    for (i in seq_len(replicates)) {
      rmatrices[[i]] <- switch(nullmodel,
        NULL_1    = NULL_1(n_rows, n_cols, n_edges),
        NULL_2    = NULL_2(n_rows, n_cols, p_rows, p_cols),
        NULL_3ROW = NULL_3ROW(n_rows, n_cols, p_rows),
        NULL_3COL = NULL_3COL(n_rows, n_cols, p_cols)
      )
    }
    rmatrices
  }

  list(NULL_1 = NULL_1, NULL_2 = NULL_2, NULL_3ROW = NULL_3ROW, NULL_3COL = NULL_3COL,
       MODEL_NAME = MODEL_NAME, NULL_MODEL = NULL_MODEL)
})
