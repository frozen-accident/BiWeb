Reading <- local({

  # File format: first line holds the quoted column labels; every following
  # line holds a quoted row label followed by that row's values.
  CREATE_FROM_MATRIX_WITH_LABELS <- function(filename) {

    if (!file.exists(filename)) {
      stop("The file could not be open. Check the name of the file.")
    }

    # With one fewer header field than data fields, read.table uses the
    # first column as row names.
    df <- read.table(filename, header = TRUE, check.names = FALSE, quote = "\"")
    matrix <- as.matrix(df)

    bipartite_web <- Bipartite(matrix)
    bipartite_web$name <- tools::file_path_sans_ext(basename(filename))
    bipartite_web
  }

  list(CREATE_FROM_MATRIX_WITH_LABELS = CREATE_FROM_MATRIX_WITH_LABELS)
})
