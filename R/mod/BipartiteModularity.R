BipartiteModularity <- setRefClass("BipartiteModularity",

  fields = list(
    matrix      = "ANY",  # Bipartite adjacency matrix
    rr          = "ANY",  # Red Nodes (rows) Communities matrix. Size = n_rows*CommunityQuantity
    tt          = "ANY",  # Blue Nodes (columns) Communities. Size = n_cols*CommunityQuantity
    n_rows      = "ANY",  # Number of rows
    n_cols      = "ANY",  # Number of columns
    n_edges     = "ANY",
    bb          = "ANY",  # Original - Null.
    index_rows  = "ANY",  # Register of the swaps in Rows.
    index_cols  = "ANY",  # Register of the swaps in Cols.
    trials      = "ANY",
    Qb          = "ANY",
    Qr          = "ANY",
    N           = "ANY",
    row_modules = "ANY",
    col_modules = "ANY",
    done        = "ANY"
  ),

  methods = list(

    initialize = function(...) {
      matrix <<- NULL; rr <<- NULL; tt <<- NULL
      n_rows <<- 0; n_cols <<- 0; n_edges <<- 0
      bb <<- NULL; index_rows <<- NULL; index_cols <<- NULL
      trials <<- 10
      Qb <<- 0; Qr <<- 0; N <<- 0
      row_modules <<- NULL; col_modules <<- NULL
      done <<- FALSE
      callSuper(...)
    },

    # Abstract
    Detect = function(...) {
      stop("Detect() is abstract; use a subclass such as LPBrim")
    },

    SortModules = function() {

      if (n_rows > n_cols) {
        idx <- order(-colSums(rr))
      } else {
        idx <- order(-colSums(tt))
      }

      index_rows <<- seq_len(n_rows)
      index_cols <<- seq_len(n_cols)

      rr <<- rr[, idx, drop = FALSE]
      tt <<- tt[, idx, drop = FALSE]

      sorted_matrix <- matrix

      row_global <- c(); col_global <- c()
      for (i in seq_len(N)) {

        row_i <- which(rr[, i] != 0)
        col_i <- which(tt[, i] != 0)

        sub <- sorted_matrix[row_i, col_i, drop = FALSE]
        row_loc <- order(-rowSums(sub))
        col_loc <- order(colSums(sub))

        row_global <- c(row_global, row_i[row_loc])
        col_global <- c(col_global, col_i[col_loc])
      }

      col_global <- rev(col_global)
      rr <<- rr[row_global, , drop = FALSE]
      tt <<- tt[col_global, , drop = FALSE]

      index_rows <<- row_global
      index_cols <<- col_global

      row_modules <<- row_modules[row_global]
      col_modules <<- col_modules[col_global]

      invisible(.self)
    }
  )
)
