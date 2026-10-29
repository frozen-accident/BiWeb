# One-sample t-test of x against mean m, mirroring MATLAB's [h p ci] = ttest(x, m).
ttest <- function(x, m) {
  if (sd(x) == 0) {
    # t.test() refuses constant data; MATLAB returns a degenerate CI instead.
    return(list(p = if (mean(x) == m) NaN else 0, ci = rep(mean(x), 2)))
  }
  r <- t.test(x, mu = m)
  list(p = r$p.value, ci = as.numeric(r$conf.int))
}

Test <- setRefClass("Test",

  fields = list(
    bipweb            = "ANY",
    nulls             = "list",
    nestvals          = "ANY",
    nestvals_rows     = "ANY",
    nestvals_cols     = "ANY",
    qb_vals           = "ANY",
    qr_vals           = "ANY",
    nest_row_contrib  = "ANY",
    nest_col_contrib  = "ANY",
    model             = "ANY",
    replicates        = "ANY",
    modul_done        = "ANY",
    nest_done         = "ANY",
    nest_contrib_done = "ANY"
  ),

  # CONSTRUCTOR AND MAIN PROCEDURE ALGORITHM
  methods = list(

    initialize = function(webbip, ...) {
      callSuper(...)
      bipweb <<- if (missing(webbip)) NULL else webbip
      nulls <<- list()
      nestvals <<- NULL; nestvals_rows <<- NULL; nestvals_cols <<- NULL
      qb_vals <<- NULL; qr_vals <<- NULL
      nest_row_contrib <<- NULL; nest_col_contrib <<- NULL
      model <<- NullModels$NULL_1
      replicates <<- 100
      modul_done <<- FALSE
      nest_done <<- FALSE
      nest_contrib_done <<- FALSE
    },

    DoNulls = function(nullmodel = NullModels$NULL_1, replic = 100) {

      modul_done <<- FALSE
      nest_done <<- FALSE

      model <<- nullmodel
      replicates <<- replic

      nulls <<- NullModels$NULL_MODEL(bipweb$adjacency, model, replicates)
      invisible(.self)
    },

    Nestedness = function() {

      if (length(nulls) == 0) {
        DoNulls()
      }
      res <- GET_DEV_NEST(bipweb, nulls)
      nestvals <<- res$out
      nestvals_rows <<- res$out_row
      nestvals_cols <<- res$out_col
      nest_done <<- TRUE
      invisible(.self)
    },

    Modularity = function() {

      if (length(nulls) == 0) {
        DoNulls()
      }
      if (!bipweb$modules$done) {
        bipweb$modules$Detect(100)
      }

      res <- GET_DEV_MODUL(bipweb, nulls)
      qb_vals <<- res$out_b
      qr_vals <<- res$out_r
      modul_done <<- TRUE
      invisible(.self)
    },

    NestednessContributions = function() {

      if (length(nulls) == 0) {
        DoNulls()
      }
      res <- GET_NEST_CONTRIBUTIONS(bipweb$adjacency, nulls)
      nest_row_contrib <<- res$c_rows
      nest_col_contrib <<- res$c_cols
      nest_contrib_done <<- TRUE
      invisible(.self)
    },

    # The methods below are static in the MATLAB version: they only use their arguments.

    GET_DEV_MODUL = function(webbip, rmatrices) {

      wQr <- webbip$modules$Qr
      wQb <- webbip$modules$Qb
      n <- length(rmatrices)

      Qb_random <- numeric(n)
      Qr_random <- numeric(n)

      for (i in seq_len(n)) {
        lp <- LPBrim(rmatrices[[i]])
        lp$Detect(10)
        Qb_random[i] <- lp$Qb
        Qr_random[i] <- lp$Qr
      }

      tb <- ttest(Qb_random, wQb)
      tr <- ttest(Qr_random, wQr)

      cib <- c(mean(tb$ci), tb$ci)
      cir <- c(mean(tr$ci), tr$ci)

      z_qb <- (wQb - mean(Qb_random)) / sd(Qb_random)
      z_qr <- (wQr - mean(Qr_random)) / sd(Qr_random)

      list(out_b = list(Qb = wQb, p = tb$p, ci = cib, zscore = z_qb),
           out_r = list(Qr = wQr, p = tr$p, ci = cir, zscore = z_qr))
    },

    GET_DEV_NEST = function(webbip, rmatrices) {

      n <- length(rmatrices)
      expect <- numeric(n)
      expect_row <- numeric(n)
      expect_col <- numeric(n)

      for (i in seq_len(n)) {
        Nodf <- NODF(rmatrices[[i]], webbip$nodf_strict)
        expect[i] <- Nodf$nodf
        expect_row[i] <- Nodf$nodf_rows
        expect_col[i] <- Nodf$nodf_cols
      }

      nest <- webbip$nestedness
      t_all <- ttest(expect, nest$nodf)
      t_row <- ttest(expect_row, nest$nodf_rows)
      t_col <- ttest(expect_col, nest$nodf_cols)

      list(
        out = list(nodf = nest$nodf, p = t_all$p, ci = c(mean(t_all$ci), t_all$ci),
                   zscore = (nest$nodf - mean(expect)) / sd(expect)),
        out_row = list(nodf = nest$nodf_rows, p = t_row$p, ci = c(mean(t_row$ci), t_row$ci),
                       zscore = (nest$nodf_rows - mean(expect_row)) / sd(expect_row)),
        out_col = list(nodf = nest$nodf_cols, p = t_col$p, ci = c(mean(t_col$ci), t_col$ci),
                       zscore = (nest$nodf_cols - mean(expect_col)) / sd(expect_col))
      )
    },

    GET_NEST_CONTRIBUTIONS = function(matrix, rmatrices) {

      N <- NODF(matrix)$nodf
      n_rows <- nrow(matrix)
      n_cols <- ncol(matrix)
      nn <- length(rmatrices)
      row_contrib <- matrix(0, nn, n_rows)
      col_contrib <- matrix(0, nn, n_cols)

      for (i in seq_len(n_rows)) {
        temp_row <- matrix[i, ]
        for (j in seq_len(nn)) {
          matrix[i, ] <- rmatrices[[j]][i, ]
          row_contrib[j, i] <- NODF(matrix, 0)$nodf
        }
        matrix[i, ] <- temp_row
      }

      for (i in seq_len(n_cols)) {
        temp_col <- matrix[, i]
        for (j in seq_len(nn)) {
          matrix[, i] <- rmatrices[[j]][, i]
          col_contrib[j, i] <- NODF(matrix, 0)$nodf
        }
        matrix[, i] <- temp_col
      }
      std_rows <- apply(row_contrib, 2, sd)
      std_cols <- apply(col_contrib, 2, sd)
      mean_rows <- colMeans(row_contrib)
      mean_cols <- colMeans(col_contrib)

      list(c_rows = (N - mean_rows) / std_rows,
           c_cols = (N - mean_cols) / std_cols)
    },

    show = function() {
      cat("Test object (null model:", NullModels$MODEL_NAME(model),
          "| replicates:", replicates, "| nulls generated:", length(nulls), ")\n")
    }
  )
)
