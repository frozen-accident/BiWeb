SpeFunc <- local({

  SP_PDI <- function(vec) {
    vec <- sort(vec, decreasing = TRUE)
    ma <- vec[1]
    R <- length(vec)
    if (R < 2) return(0)
    sum((vec[1] / ma - vec[-1] / ma) / (R - 1))
  }

  SPECIFICITY <- function(webmatrix) {
    apply(webmatrix, 1, SP_PDI)
  }

  RESOURCE_RANGE <- function(webmatrix) {
    row_degrees <- rowSums(webmatrix > 0)
    R <- ncol(webmatrix)
    (R - row_degrees) / (R - 1)
  }

  SPECIES_SPECIFICITY_INDEX <- function(webmatrix) {
    R <- ncol(webmatrix)
    norm_factor <- R * sqrt((R - 1) / R)
    Pbar <- rowMeans(webmatrix)
    tspe <- rowSums((webmatrix - Pbar)^2)
    sqrt(tspe) / (norm_factor * Pbar)
  }

  IR <- function(webmatrix) {
    n_rows <- nrow(webmatrix)
    # Correlation matrix of responses
    rho <- suppressWarnings(cor(t(webmatrix)))
    # Standard deviations of responses (normalized by N, as std(x,1))
    sig <- sqrt(rowMeans((webmatrix - rowMeans(webmatrix))^2))
    # Correction factor
    cfac <- n_rows * (n_rows - 1)
    # Measures
    tR <- 0
    tI <- 0
    for (i in seq_len(n_rows - 1)) {
      for (j in colon(i + 1, n_rows)) {
        tR <- tR + (sig[i] - sig[j])^2
        tI <- tI + sig[i] * sig[j] * (1 - rho[i, j])
      }
    }

    R <- tR / (2 * cfac)
    I <- tI / cfac

    VarW <- mean((webmatrix - mean(webmatrix))^2)

    if (VarW == 0) {
      R <- 0
      I <- 0
    } else {
      R <- R / VarW
      I <- I / VarW
    }

    list(R = R, I = I)
  }

  list(SP_PDI = SP_PDI, SPECIFICITY = SPECIFICITY, RESOURCE_RANGE = RESOURCE_RANGE,
       SPECIES_SPECIFICITY_INDEX = SPECIES_SPECIFICITY_INDEX, IR = IR)
})
