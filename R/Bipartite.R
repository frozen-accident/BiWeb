#==========================================================================
# Name        : Bipartite.R
# Author      : Cesar Flores (MATLAB original)
# Description : Represent a complete data structure of a bipartite network.
#==========================================================================

Bipartite <- setRefClass("Bipartite",

  fields = list(
    webmatrix      = "ANY", # Interaction matrix (not necessarily a binary matrix).
    adjacency      = "ANY", # Adjacency matrix (binary matrix)
    modules        = "ANY",
    nestedness     = "ANY",
    num_edges      = "ANY", # Number of edges
    connectance    = "ANY", # Fill of the webmatrix
    n_rows         = "ANY", # Number of rows
    n_cols         = "ANY", # Number of columns
    size_webmatrix = "ANY",
    specificity    = "ANY",
    row_degrees    = "ANY",
    col_degrees    = "ANY",
    rr             = "ANY",
    ssi            = "ANY",
    bperf          = "ANY",
    name           = "ANY",
    tests          = "ANY",
    nodf_strict    = "ANY",
    row_labels     = "ANY",
    col_labels     = "ANY",
    printer        = "ANY",
    plotter        = "ANY"
  ),

  methods = list(

    # web: a numeric/logical matrix (row/column names are used as labels if
    # present) or the path of a whitespace-delimited text file in matrix format.
    initialize = function(web, namebip) {

      if (missing(web)) {
        stop("You need to specify a double matrix or a txt file in matrix format")
      }

      if (is.character(web)) {
        if (missing(namebip)) namebip <- tools::file_path_sans_ext(basename(web))
        web <- as.matrix(read.table(web))
        dimnames(web) <- NULL
      } else if (missing(namebip)) {
        namebip <- "No name"
      }

      web <- 1.0 * as.matrix(web)
      keep_rows <- rowSums(web) != 0
      keep_cols <- colSums(web) != 0

      web <- web[keep_rows, keep_cols, drop = FALSE]
      labels_r <- rownames(web)
      labels_c <- colnames(web)
      dimnames(web) <- NULL

      name <<- namebip
      webmatrix <<- web

      # General Info
      n_rows <<- nrow(web)
      n_cols <<- ncol(web)
      size_webmatrix <<- n_rows * n_cols

      # Connectance
      adjacency <<- webmatrix > 0
      num_edges <<- sum(adjacency)
      connectance <<- sum(adjacency) / length(adjacency)

      # Specificity and all
      row_degrees <<- rowSums(adjacency)
      col_degrees <<- colSums(adjacency)
      specificity <<- SpeFunc$SPECIFICITY(webmatrix)
      rr <<- SpeFunc$RESOURCE_RANGE(webmatrix)
      ssi <<- SpeFunc$SPECIES_SPECIFICITY_INDEX(webmatrix)
      bperf <<- apply(webmatrix, 2, max)

      nodf_strict <<- 1
      nestedness <<- NODF(adjacency, nodf_strict)
      modules <<- LPBrim(webmatrix)
      tests <<- Test(.self)
      printer <<- Printer(.self)
      plotter <<- PlotWebs(.self)

      row_labels <<- if (is.null(labels_r)) sprintf("row%03i", seq_len(n_rows)) else labels_r
      col_labels <<- if (is.null(labels_c)) sprintf("col%03i", seq_len(n_cols)) else labels_c
    },

    # The default show() would recurse through tests/printer/plotter, which point back here.
    show = function() {
      cat(sprintf("Bipartite network '%s': %d rows x %d cols, %d interactions\n",
                  name, n_rows, n_cols, as.integer(num_edges)))
      cat(sprintf("  connectance %.3f | NODF %.3f (rows %.3f, cols %.3f)\n",
                  connectance, nestedness$nodf, nestedness$nodf_rows, nestedness$nodf_cols))
      if (modules$done) {
        cat(sprintf("  modules %d | Qb %.3f | Qr %.3f\n", modules$N, modules$Qb, modules$Qr))
      }
    }
  )
)
