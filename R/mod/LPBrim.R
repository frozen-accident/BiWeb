#==========================================================================
# Name        : LPBrim.R
# Author      : Cesar Flores (MATLAB original)
# Description : Represents the LP&BRIM algorithm
#==========================================================================

LPBrim <- setRefClass("LPBrim",

  contains = "BipartiteModularity",

  fields = list(
    red_labels  = "ANY",
    blue_labels = "ANY"
  ),

  # CONSTRUCTOR AND MAIN PROCEDURES ALGORITHM
  methods = list(

    initialize = function(bipmatrix, ...) {

      callSuper(...)
      red_labels <<- 0
      blue_labels <<- 0
      if (missing(bipmatrix)) return(invisible(.self))

      matrix <<- unname(as.matrix(bipmatrix)) > 0

      n_rows <<- nrow(matrix)
      n_cols <<- ncol(matrix)

      trials <<- 100

      if (all(!matrix)) {
        bb <<- 1.0 * matrix
      } else {
        coldeg <- colSums(matrix)
        rowdeg <- rowSums(matrix)
        n_edges <<- sum(rowdeg)
        bb <<- matrix - (1 / n_edges) * outer(rowdeg, coldeg)
      }
    },

    Detect = function(ntrials) {

      if (!missing(ntrials)) {
        trials <<- ntrials
      }

      LP()

      BRIM()

      CleanCommunities()

      CalculateQrValue()

      idx <- which(rr != 0, arr.ind = TRUE)
      row_modules <<- unname(idx[order(idx[, 1]), 2])

      idx <- which(tt != 0, arr.ind = TRUE)
      col_modules <<- unname(idx[order(idx[, 1]), 2])

      SortModules()

      done <<- TRUE

      invisible(.self)
    },

    LP = function() {

      red_labels <<- seq_len(n_rows)
      blue_labels <<- seq_len(n_cols)
      cmax <- max(n_rows, n_cols)
      rr <<- diag(1, n_rows, cmax)
      tt <<- diag(1, n_cols, cmax)

      n_edges <<- sum(matrix)

      qmax <- CalculateQValue()
      qmaxglobal <- qmax
      ttMaxGlobal <- tt
      rrMaxGlobal <- rr
      for (trial in seq_len(trials)) {

        cmax <- max(n_rows, n_cols)
        red_labels <<- seq_len(n_rows)
        blue_labels <<- seq_len(n_cols)
        rr <<- diag(1, n_rows, cmax)
        tt <<- diag(1, n_cols, cmax)
        qmax <- CalculateQValue()
        rrMax <- rr
        ttMax <- tt

        repeat {

          # Propagate red labels to blue labels
          for (i in seq_len(n_cols)) {
            labels <- red_labels[matrix[, i]]
            if (length(labels) > 0) {
              uniq <- sort(unique(labels))
              uniq <- uniq[sample.int(length(uniq))]
              count <- vapply(uniq, function(x) sum(x == labels), numeric(1))
              blue_labels[i] <<- uniq[which.max(count)]
            }
          }

          uniq <- sort(unique(blue_labels))
          nbluec <- length(uniq)
          blue_labels <<- match(blue_labels, uniq)

          # Propagate blue labels to red labels
          for (i in seq_len(n_rows)) {
            labels <- blue_labels[matrix[i, ]]
            if (length(labels) > 0) {
              uniq <- sort(unique(labels))
              uniq <- uniq[sample.int(length(uniq))]
              count <- vapply(uniq, function(x) sum(x == labels), numeric(1))
              red_labels[i] <<- uniq[which.max(count)]
            }
          }

          uniq <- sort(unique(red_labels))
          nredc <- length(uniq)
          red_labels <<- match(red_labels, uniq)

          cmax <- max(nbluec, nredc)

          rr <<- matrix(0, n_rows, cmax)
          tt <<- matrix(0, n_cols, cmax)
          rr[cbind(seq_len(n_rows), red_labels)] <<- 1
          tt[cbind(seq_len(n_cols), blue_labels)] <<- 1

          q <- CalculateQValue()

          if (isTRUE(q > qmax)) {
            qmax <- q
            rrMax <- rr
            ttMax <- tt
          } else {
            break
          }
        }

        if (isTRUE(qmax > qmaxglobal)) {
          qmaxglobal <- qmax
          rrMaxGlobal <- rrMax
          ttMaxGlobal <- ttMax
        }
      }
      Qb <<- qmaxglobal
      rr <<- rrMaxGlobal
      tt <<- ttMaxGlobal
      invisible(.self)
    },

    BRIM = function() {

      inducingBlueFlag <- TRUE
      n_edges <<- sum(matrix)

      preQ <- CalculateQValue()
      repeat {

        if (inducingBlueFlag) {
          AssignBlueNodes()
        } else {
          AssignRedNodes()
        }

        newQ <- CalculateQValue()

        if (!isTRUE(newQ > preQ)) {
          break
        } else {
          preQ <- newQ
          inducingBlueFlag <- !inducingBlueFlag
        }
      }
      invisible(.self)
    },

    CalculateQValue = function() {
      sum(diag(t(rr) %*% bb %*% tt)) / n_edges
    },

    CalculateQrValue = function() {

      Qr <<- 0

      for (i in seq_len(N)) {
        row_index <- which(rr[, i] != 0)
        col_index <- which(tt[, i] != 0)
        Qr <<- Qr + sum(matrix[row_index, col_index, drop = FALSE] > 0)
      }

      Qr <<- Qr / n_edges
      invisible(.self)
    },

    AssignRedNodes = function() {
      maxind <- max.col(bb %*% tt, ties.method = "first")
      modules <- diag(ncol(tt))
      rr <<- modules[maxind, , drop = FALSE]
      invisible(.self)
    },

    AssignBlueNodes = function() {
      maxind <- max.col(t(bb) %*% rr, ties.method = "first")
      modules <- diag(ncol(rr))
      tt <<- modules[maxind, , drop = FALSE]
      invisible(.self)
    },

    CleanCommunities = function() {

      # CLEAN THE COMMUNITIES
      com1 <- which(colSums(tt != 0) > 0)
      com2 <- which(colSums(rr != 0) > 0)

      if (length(com1) >= length(com2)) { # Not sure if this validation is necessary.
        rr <<- rr[, com1, drop = FALSE]
        tt <<- tt[, com1, drop = FALSE]
      } else {
        rr <<- rr[, com2, drop = FALSE]
        tt <<- tt[, com2, drop = FALSE]
      }

      N <<- ncol(rr)
      invisible(.self)
    }
  )
)
