Printer <- setRefClass("Printer",

  fields = list(
    bipweb = "ANY"
  ),

  methods = list(

    initialize = function(webbip, ...) {
      callSuper(...)
      bipweb <<- if (missing(webbip)) NULL else webbip
    },

    SpeciesLevel = function(toscreen = TRUE, tofile = TRUE) {

      webbip <- bipweb
      namebip <- gsub(" ", "_", webbip$name)
      filename <- paste0(namebip, "-sp.txt")

      colwidth <- 9
      formatstring <- paste0("%-", colwidth, "s")
      formatint <- paste0("%-", colwidth, "d")
      formatfloat <- paste0("%-", colwidth, ".3f")
      estring <- "----"

      header <- paste0(sprintf(formatstring, c("Specie", "Level", "Degree", "SPE", "SSI", "RR")), collapse = "")
      if (webbip$modules$done) header <- paste0(header, sprintf(formatstring, "Module"))
      header <- paste0(header, "\n")

      row_names <- substr(gsub(" ", "_", webbip$row_labels), 1, colwidth)
      rows <- paste0(sprintf(formatstring, row_names),
                     sprintf(formatstring, "top"),
                     sprintf(formatint, as.integer(webbip$row_degrees)),
                     sprintf(formatfloat, webbip$specificity),
                     sprintf(formatfloat, webbip$ssi),
                     sprintf(formatfloat, webbip$rr))
      if (webbip$modules$done) rows <- paste0(rows, sprintf(formatint, as.integer(webbip$modules$row_modules)))

      col_names <- substr(gsub(" ", "_", webbip$col_labels), 1, colwidth)
      cols <- paste0(sprintf(formatstring, col_names),
                     sprintf(formatstring, "bottom"),
                     sprintf(formatint, as.integer(webbip$col_degrees)),
                     strrep(sprintf(formatstring, estring), 3))
      if (webbip$modules$done) cols <- paste0(cols, sprintf(formatint, as.integer(webbip$modules$col_modules)))

      info <- paste0(c(rows, cols), "\n", collapse = "")

      if (toscreen) {
        cat(header, info, sep = "")
      }

      if (tofile) {
        cat(header, info, sep = "", file = filename)
      }

      invisible(.self)
    },

    NetworkLevelSingle = function(toscreen = TRUE, tofile = TRUE,
                                  filename = paste0(bipweb$name, "_networklevel.txt")) {
      # To complete, maybe have to delete at the end of the day
      webbip <- bipweb

      fmt <- function(x) sprintf("%5.2f", x)
      fmtint <- function(x) sprintf("%5d", as.integer(x))
      header <- ""
      if (!file.exists(filename)) {
        header <- paste0(header, "\nGeneral properties\n")
        header <- paste0(header, "------------------\n")
        header <- paste0(header, "\tNetwork Name:      ", webbip$name, "\n")
        header <- paste0(header, "\tRows:              ", fmtint(webbip$n_rows), "\n")
        header <- paste0(header, "\tCols:              ", fmtint(webbip$n_cols), "\n")
        header <- paste0(header, "\tSize:              ", fmtint(webbip$size_webmatrix), "\n")
        header <- paste0(header, "\tInteractions:      ", fmtint(webbip$num_edges), "\n")
        header <- paste0(header, "\tConnectance:       ", fmt(webbip$connectance), "\n")
        header <- paste0(header, "\t<Specificity>:     ", fmt(mean(webbip$specificity)), "\n")
        header <- paste0(header, "\t<Res. Range.>:     ", fmt(mean(webbip$rr)), "\n")
        header <- paste0(header, "\tNODF:              ", fmt(webbip$nestedness$nodf), "\n")
        header <- paste0(header, "\tNODF Rows:         ", fmt(webbip$nestedness$nodf_rows), "\n")
        header <- paste0(header, "\tNODF Cols:         ", fmt(webbip$nestedness$nodf_cols), "\n")
      }

      if (toscreen) {
        cat(header)
      }

      invisible(.self)
    },

    NetworkLevel = function(toscreen = TRUE, tofile = TRUE) {
      webbip <- bipweb
      tests <- webbip$tests
      modules <- webbip$modules
      filename <- "networklevel.txt"

      colwidth <- 9
      formatstring <- paste0("%-", colwidth, "s")
      formatint <- paste0("%-", colwidth, "d")
      formatfloat <- paste0("%-", colwidth, ".3f")
      estring <- "----"
      fs <- function(x) paste0(sprintf(formatstring, x), collapse = "")
      ff <- function(x) paste0(sprintf(formatfloat, x), collapse = "")
      fi <- function(x) sprintf(formatint, as.integer(x))

      header <- fs(c("Name", "Conn", "Size", "n_rows", "n_cols",
                     "nodf", "nodf_r", "nodf_c", "aspe", "arr",
                     "resp", "inc", "mN", "mQb", "mQr",
                     "n_sim", "n_pval", "nicLow", "nicUp", "n_z", "nco_sim",
                     "nco_pval", "ncoicLow", "ncoicUp", "nco_z", "nro_sim",
                     "nro_pval", "nroicLow", "nroicUp", "nro_z",
                     "mQr_sim", "mQr_val", "mQrIclow", "mQrIcup", "mQr_z",
                     "mQb_sim", "mQb_val", "mQbIclow", "mQbIcup", "mQb_z",
                     "tReps", "null_name"))
      header <- paste0(header, "\n")

      ir <- SpeFunc$IR(webbip$webmatrix)
      namebip <- gsub(" ", "_", webbip$name)
      info <- sprintf(formatstring, substr(namebip, 1, colwidth - 1))
      info <- paste0(info, ff(webbip$connectance),
                     fi(webbip$size_webmatrix), fi(webbip$n_rows), fi(webbip$n_cols),
                     ff(c(webbip$nestedness$nodf, webbip$nestedness$nodf_rows, webbip$nestedness$nodf_cols,
                          mean(webbip$specificity), mean(webbip$rr), ir$R, ir$I)))

      if (modules$done) {
        info <- paste0(info, fi(modules$N), ff(c(modules$Qb, modules$Qr)))
      } else {
        info <- paste0(info, fs(rep(estring, 3)))
      }

      # simulated mean, p-value, CI low, CI up, z-score
      stat_cols <- function(v) ff(c(v$ci[1], v$p, v$ci[2], v$ci[3], v$zscore))

      if (tests$nest_done) {
        info <- paste0(info, stat_cols(tests$nestvals),
                       stat_cols(tests$nestvals_cols),
                       stat_cols(tests$nestvals_rows))
      } else {
        info <- paste0(info, fs(rep(estring, 12)))
      }

      if (tests$modul_done) {
        info <- paste0(info, stat_cols(tests$qr_vals), stat_cols(tests$qb_vals))
      } else {
        info <- paste0(info, fs(rep(estring, 8)))
      }

      if (tests$modul_done || tests$nest_done) {
        info <- paste0(info, fi(tests$replicates),
                       sprintf(formatstring, paste0("NullModels.", NullModels$MODEL_NAME(tests$model))))
      } else {
        info <- paste0(info, fs(rep(estring, 2)))
      }

      info <- paste0(info, "\n")

      if (toscreen) {
        cat(header, info, sep = "")
      }

      if (tofile) {
        if (!file.exists(filename)) {
          cat(header, info, sep = "", file = filename)
        } else {
          cat(info, file = filename, append = TRUE)
        }
      }

      invisible(.self)
    },

    show = function() {
      cat("Printer for network '", bipweb$name, "'\n", sep = "")
    }
  )
)
