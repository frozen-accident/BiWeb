MatrixNull <- local({

  NESTED_TO_MODULAR <- function(n_size, n_links) {
    n <- n_size
    matrix <- matrix(0, n, n)
    matrix[lower.tri(matrix, diag = TRUE)] <- 1
    matrix <- matrix[n:1, , drop = FALSE]
    nested_sites <- c()
    modular_sites <- c()
    for (i in seq_len(floor(n / 2))) {
      nested_sites <- c(nested_sites, colon((i - 1) * n + n / 2 + 1, n * i - i + 1))
      nested_sites <- c(nested_sites, colon(n / 2 * n + (i - 1) * n + 1, n / 2 * n + (i - 1) * n + 1 + n / 2 - i))
      modular_sites <- c(modular_sites, colon(n / 2 * n + (i - 1) * n + n / 2 + 1, n / 2 * n + (i - 1) * n + n))
    }

    n_nest <- length(nested_sites)
    n_modul <- length(modular_sites)

    random_nest_switched <- nested_sites[sample.int(n_nest, min(n_nest, n_links))]
    random_modul_switched <- modular_sites[sample.int(n_modul, min(n_modul, n_links))]

    matrix[random_nest_switched] <- 0
    matrix[random_modul_switched] <- 1

    imagesc(matrix)

    list(matrix = matrix, nested_sites = nested_sites, modular_sites = modular_sites)
  }

  SORT_MATRIX <- function(matrix) {
    ir <- order(-rowSums(matrix))
    ic <- order(-colSums(matrix))
    matrix[ir, ic, drop = FALSE]
  }

  RANDOM_SORT <- function(matrix) {
    matrix[sample.int(nrow(matrix)), sample.int(ncol(matrix)), drop = FALSE]
  }

  list(NESTED_TO_MODULAR = NESTED_TO_MODULAR, SORT_MATRIX = SORT_MATRIX,
       RANDOM_SORT = RANDOM_SORT)
})
