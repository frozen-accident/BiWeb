PlotWebs <- setRefClass("PlotWebs",

  fields = list(
    pos_x                 = "ANY", # X coordinate of each cell in a unit square
    pos_y                 = "ANY", # Y coordinate of each cell in a unit square
    cell_color            = "ANY", # Cell color
    back_color            = "ANY", # Back color
    margin                = "ANY", # Margin between cells
    ax                    = "ANY", # Axis ("image" keeps a 1:1 aspect ratio)
    border_color          = "ANY",
    border_line_width     = "ANY",

    matrix                = "ANY",
    biweb                 = "ANY",
    row_labels            = "ANY",
    col_labels            = "ANY",
    use_labels            = "ANY",
    FontSize              = "ANY", # In points
    index_rows            = "ANY",
    index_cols            = "ANY",
    colors                = "ANY",
    n_rows                = "ANY", # row number
    n_cols                = "ANY", # col number

    radius                = "ANY",
    vertical_margin       = "ANY",
    horizontal_proportion = "ANY",
    row_pos               = "ANY",
    col_pos               = "ANY",
    bead_color            = "ANY",
    link_color            = "ANY"
  ),

  methods = list(

    initialize = function(bipmatrix_or_biweb, ...) {

      callSuper(...)
      cell_color <<- "white"
      back_color <<- rgb(0, 0, 128, maxColorValue = 255)
      margin <<- 0.12
      ax <<- "image"
      border_color <<- NA
      border_line_width <<- 0.0005

      biweb <<- NULL
      row_labels <<- character(0)
      col_labels <<- character(0)
      use_labels <<- TRUE
      FontSize <<- 5
      index_rows <<- NULL
      index_cols <<- NULL

      radius <<- 0.5
      vertical_margin <<- 0.12
      horizontal_proportion <<- 0.5
      bead_color <<- "black"
      link_color <<- "black"

      if (missing(bipmatrix_or_biweb)) return(invisible(.self))

      if (inherits(bipmatrix_or_biweb, "Bipartite")) {
        biweb <<- bipmatrix_or_biweb
        matrix <<- biweb$adjacency
      } else {
        matrix <<- bipmatrix_or_biweb
      }

      n_rows <<- nrow(matrix)
      n_cols <<- ncol(matrix)

      colors <<- jet(64)[c(23,2,13,42,57,20,15,11,9,16,3,28,26,24,46,59,41,18,56,40,17,48,27,53,6,62,5,60,14,32,64,19,
                           36,58,39,21,4,8,35,30,50,63,25,51,55,34,61,37,47,44,54,43,38,12,52,33,31,1,22,29,10,45,49,7)]

      FillPositions()
    },

    FillPositions = function() {
      pos_x <<- matrix(seq_len(n_cols), n_rows, n_cols, byrow = TRUE)
      pos_y <<- matrix((n_rows + 1) - seq_len(n_rows), n_rows, n_cols)

      maxd <- max(n_rows, n_cols)
      row_pos <<- seq(maxd, 1, length.out = n_rows)
      col_pos <<- seq(maxd, 1, length.out = n_cols)
      invisible(.self)
    },

    # Color assigned to each module number (cycles through the palette).
    ModuleColor = function(m) {
      colors[(m - 1) %% length(colors) + 1]
    },

    PlotAspect = function() {
      if (ax %in% c("image", "equal")) 1 else NA
    },

    SetMargins = function() {
      par(mar = if (use_labels) c(3, 3, 1, 1) else c(1, 1, 1, 1))
    },

    PlotBMatrix = function() {
      op <- SetMargins(); on.exit(par(op))
      maxd <- max(n_rows, n_cols)
      x1 <- 1; x2 <- 1 + maxd * horizontal_proportion

      index_rows <<- seq_len(n_rows)
      index_cols <<- seq_len(n_cols)

      local_matrix <- matrix

      StartBPlot()
      links <- which(local_matrix == 1, arr.ind = TRUE)
      segments(1, row_pos[links[, 1]], x2, col_pos[links[, 2]], col = "black")

      DrawCircle(rep(x1, n_rows), row_pos, "black")
      DrawCircle(rep(x2, n_cols), col_pos, "black")

      ApplyBasicBFormat()
    },

    PlotBNestedMatrix = function() {
      op <- SetMargins(); on.exit(par(op))
      maxd <- max(n_rows, n_cols)
      x1 <- 1; x2 <- 1 + maxd * horizontal_proportion

      index_rows <<- order(-rowSums(matrix))
      index_cols <<- order(-colSums(matrix))

      local_matrix <- matrix[index_rows, index_cols, drop = FALSE]

      StartBPlot()
      links <- which(local_matrix == 1, arr.ind = TRUE)
      segments(1, row_pos[links[, 1]], x2, col_pos[links[, 2]], col = link_color)

      DrawCircle(rep(x1, n_rows), row_pos, bead_color)
      DrawCircle(rep(x2, n_cols), col_pos, bead_color)

      ApplyBasicBFormat()
    },

    PlotBModularMatrix = function() {
      op <- SetMargins(); on.exit(par(op))
      maxd <- max(n_rows, n_cols)
      x1 <- 1; x2 <- 1 + maxd * horizontal_proportion

      if (is.null(biweb)) {
        biweb <<- Bipartite(matrix)
      }
      if (!biweb$modules$done) {
        biweb$modules$Detect()
      }

      index_rows <<- biweb$modules$index_rows
      index_cols <<- rev(biweb$modules$index_cols)

      row_mod <- biweb$modules$row_modules
      col_mod <- rev(biweb$modules$col_modules)

      local_matrix <- matrix[index_rows, index_cols, drop = FALSE]

      StartBPlot()
      links <- which(local_matrix == 1, arr.ind = TRUE)
      same <- row_mod[links[, 1]] == col_mod[links[, 2]]
      link_cols <- ifelse(same, ModuleColor(row_mod[links[, 1]]), "black")
      segments(1, row_pos[links[, 1]], x2, col_pos[links[, 2]], col = link_cols)

      DrawCircle(rep(x1, n_rows), row_pos, ModuleColor(row_mod))
      DrawCircle(rep(x2, n_cols), col_pos, ModuleColor(col_mod))

      ApplyBasicBFormat()
    },

    DrawLine = function(start_cord, end_cord, color) {
      segments(start_cord[1], start_cord[2], end_cord[1], end_cord[2], col = color, lwd = 1.0)
    },

    DrawCircle = function(x, y, color) {
      if (length(x) == 0) return(invisible(.self))
      r <- radius - vertical_margin
      symbols(x, y, circles = rep(r, length(x)), inches = FALSE, add = TRUE,
              bg = color, fg = color)
      invisible(.self)
    },

    PlotMatrix = function() {
      op <- SetMargins(); on.exit(par(op))
      index_rows <<- seq_len(n_rows)
      index_cols <<- seq_len(n_cols)

      StartMatrixPlot()
      cells <- which(matrix > 0, arr.ind = TRUE)
      DrawCell(cells[, 1], cells[, 2], cell_color)

      ApplyBasicFormat()
    },

    PlotNestedMatrix = function() {
      op <- SetMargins(); on.exit(par(op))
      index_rows <<- order(-rowSums(matrix))
      index_cols <<- order(-colSums(matrix))

      local_matrix <- matrix[index_rows, index_cols, drop = FALSE]

      StartMatrixPlot()
      cells <- which(local_matrix > 0, arr.ind = TRUE)
      DrawCell(cells[, 1], cells[, 2], cell_color)

      ApplyBasicFormat()
    },

    PlotModularMatrix = function() {
      op <- SetMargins(); on.exit(par(op))
      if (is.null(biweb)) {
        biweb <<- Bipartite(matrix)
      }
      if (!biweb$modules$done) {
        biweb$modules$Detect()
      }

      index_rows <<- biweb$modules$index_rows
      index_cols <<- biweb$modules$index_cols

      row_mod <- biweb$modules$row_modules
      col_mod <- biweb$modules$col_modules

      local_matrix <- matrix[index_rows, index_cols, drop = FALSE]

      StartMatrixPlot()
      cells <- which(local_matrix > 0, arr.ind = TRUE)
      same <- row_mod[cells[, 1]] == col_mod[cells[, 2]]
      DrawCell(cells[, 1], cells[, 2], ifelse(same, ModuleColor(row_mod[cells[, 1]]), cell_color))

      ApplyBasicFormat()
    },

    MatrixLimits = function() {
      c(0.5 - 1.1 * margin, 0.5 + n_cols + margin, 0.5 - margin, 0.5 + n_rows + margin)
    },

    StartMatrixPlot = function() {
      lim <- MatrixLimits()
      plot.new()
      plot.window(xlim = lim[1:2], ylim = lim[3:4], asp = PlotAspect(), xaxs = "i", yaxs = "i")
      DrawBack()
    },

    ApplyBasicFormat = function() {

      if (use_labels) {

        FillLabels()

        cex <- FontSize / par("ps")
        text(0, pos_y[, 1], row_labels, adj = c(1, 0.5), cex = cex, xpd = NA)
        text(pos_x[1, ], 0, col_labels, adj = c(1, 0.5), srt = 90, cex = cex, xpd = NA)
      }

      lim <- MatrixLimits()
      rect(lim[1], lim[3], lim[2], lim[4], border = "black", xpd = NA)
      invisible(.self)
    },

    StartBPlot = function() {
      maxd <- max(n_rows, n_cols)
      x1 <- 1; x2 <- 1 + maxd * horizontal_proportion
      plot.new()
      plot.window(xlim = c(x1 - radius, x2 + radius),
                  ylim = c(1 - radius, maxd + radius),
                  asp = PlotAspect(), xaxs = "i", yaxs = "i")
    },

    ApplyBasicBFormat = function() {
      maxd <- max(n_rows, n_cols)
      x1 <- 1; x2 <- 1 + maxd * horizontal_proportion

      if (use_labels) {

        FillLabels()

        cex <- FontSize / par("ps")
        text(x1 - radius, row_pos, row_labels, adj = c(1, 0.5), cex = cex, xpd = NA)
        text(x2 + radius, col_pos, col_labels, adj = c(0, 0.5), cex = cex, xpd = NA)
      }
      invisible(.self)
    },

    FillLabels = function() {

      if (is.null(biweb)) {
        row_labels <<- sprintf("row%03i", seq_len(n_rows))
        col_labels <<- sprintf("col%03i", seq_len(n_cols))
      } else {
        row_labels <<- biweb$row_labels
        col_labels <<- biweb$col_labels
      }

      row_labels <<- row_labels[index_rows]
      col_labels <<- col_labels[index_cols]
      invisible(.self)
    },

    DrawBack = function() {
      lim <- MatrixLimits()
      rect(lim[1], lim[3], lim[2], lim[4], col = back_color, border = NA)
      invisible(.self)
    },

    DrawCell = function(i, j, color) {
      x <- pos_x[cbind(i, j)]
      y <- pos_y[cbind(i, j)]
      rect(x - 0.5 + margin, y - 0.5 + margin, x + 0.5 - margin, y + 0.5 - margin,
           col = color, border = NA)
      invisible(.self)
    },

    show = function() {
      cat(sprintf("PlotWebs for a %d x %d matrix\n", n_rows, n_cols))
    }
  )
)
