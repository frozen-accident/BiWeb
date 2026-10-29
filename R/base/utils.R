# Small helpers that stand in for MATLAB built-ins used across the toolbox.

# MATLAB's jet(m) colormap, returned as a vector of hex colors.
jet <- function(m = 64) {
  n <- ceiling(m / 4)
  u <- c((1:n) / n, rep(1, n - 1), (n:1) / n)
  g <- ceiling(n / 2) - (m %% 4 == 1) + seq_along(u)
  r <- g + n
  b <- g - n
  g <- g[g <= m]; r <- r[r <= m]; b <- b[b >= 1]
  J <- matrix(0, m, 3)
  J[r, 1] <- u[seq_along(r)]
  J[g, 2] <- u[seq_along(g)]
  J[b, 3] <- u[(length(u) - length(b) + 1):length(u)]
  rgb(J[, 1], J[, 2], J[, 3])
}

# MATLAB's imagesc(matrix): row 1 at the top, jet colormap.
imagesc <- function(m) {
  m <- 1.0 * as.matrix(m)
  nr <- nrow(m); nc <- ncol(m)
  image(1:nc, 1:nr, t(m[nr:1, , drop = FALSE]), col = jet(64),
        xlab = "", ylab = "", axes = FALSE, asp = 1)
  box()
}

# MATLAB's a:b, which is empty (rather than decreasing) when b < a.
colon <- function(a, b) if (b < a) numeric(0) else seq(a, b)
