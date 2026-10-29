# Load the code (the R equivalent of addpath(genpath('matlab/')))
# Locate the folder holding this script, so it works from any working directory
# (via source() or Rscript); otherwise fall back to the working directory.
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) dirname(sub("^--file=", "", arg)) else getwd()
})
r_files <- list.files(file.path(here, "R"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
if (length(r_files) == 0) {
  stop("Could not find the R/ folder in '", here, "'. Set the working directory to the BiWeb folder first, e.g. setwd(\"path/to/BiWeb\").")
}
for (f in r_files) source(f)

# Create an object of the class Bipartite. By just creating it,
# nestedness will be calculated;
w <- Reading$CREATE_FROM_MATRIX_WITH_LABELS(file.path(here, "sampledata", "rodents.web"))

# Accessing the values of nestedness
nodf <- w$nestedness$nodf
nodf_up <- w$nestedness$nodf_rows
nodf_low <- w$nestedness$nodf_cols

cat("\nValues of Nestedness using NODF algorithm:\n")
cat(sprintf("\tNODF:                 \t%f\n", nodf))
cat(sprintf("\tNODF in rows:         \t%f\n", nodf_up))
cat(sprintf("\tNODF in columns:      \t%f\n", nodf_low))

# Accessing the values of modularity using LP&BRIM algorithm
w$modules$Detect()
N <- w$modules$N
Qb <- w$modules$Qb
Qr <- w$modules$Qr

cat("\nValues of Modularity using LP&BRIM algorithm:\n")
cat(sprintf("\tNumber of modules:    \t%f\n", N))
cat(sprintf("\tBipartite modularity: \t%f\n", Qb))
cat(sprintf("\tRealized modularity:  \t%f\n", Qr))


# Print species level properties;
# The properties will be printed to screen and to file rodents-sp.txt
w$printer$SpeciesLevel()

# Create statistics for modularity and nestedness;
n_trials <- 10 # Increase this value for a more accurate result
w$tests$DoNulls(NullModels$NULL_1, n_trials) # Using bernoulli null model
w$tests$Nestedness()
w$tests$Modularity()
w$printer$NetworkLevel() # Add a row with the statistics to networklevel file
w$tests$DoNulls(NullModels$NULL_2, n_trials) # Using average in columns+rows
w$tests$Nestedness()
w$tests$Modularity()
w$printer$NetworkLevel()
w$tests$DoNulls(NullModels$NULL_3COL, n_trials) # Using average in columns
w$tests$Nestedness()
w$tests$Modularity()
w$printer$NetworkLevel()
w$tests$DoNulls(NullModels$NULL_3ROW, n_trials) # Using average in rows
w$tests$Nestedness()
w$tests$Modularity()
w$printer$NetworkLevel()

# Plot network as matrices
# (each call draws a new plot; when run with Rscript they become pages of Rplots.pdf)
w$plotter$PlotMatrix()
w$plotter$PlotNestedMatrix()
w$plotter$PlotModularMatrix()

# Plot network as beads (graph)
w$plotter$PlotBMatrix()
w$plotter$PlotBNestedMatrix()
w$plotter$PlotBModularMatrix()
