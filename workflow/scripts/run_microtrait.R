#!/usr/bin/Rscript

# https://github.com/ukaraoz/microtrait

## if library errors happen uncomment the below
##.libPaths(c("/home/susbus/R/x86_64-pc-linux-gnu-library/4.4",
##            "/usr/local/lib/R/site-library",
##            "/usr/lib/R/site-library",
##            "/usr/lib/R/library")) 


############################## LOG
sink(file=file(snakemake@log[["out"]], open="wt"), type=c("output", "message"))

# Set TMPDIR using R's Sys.setenv
# Sys.setenv(TMPDIR = dirname(snakemake@output[["rds"]]))

############################## LIBS
# loading library
library(microtrait)
library(tictoc)

setequal<- Biostrings::setequal
intersect<-Biostrings::intersect
union<-Biostrings::union
collapse<-Biostrings::collapse
setdiff<-Biostrings::setdiff

# running microtrait

# Debugging: Print the value of snakemake@input$mags to check the path
message("Input MAGS: ", snakemake@input[["mags"]], "\n")

# Use the path directly from snakemake input
genomes_dir = snakemake@input[["mags"]]

# Debugging: Check if the directory is correctly specified
message("Genomes directory: ", genomes_dir, "\n")

genomes_files = list.files(genomes_dir, full.names = TRUE, recursive = TRUE, pattern = ".fa$")

# Debugging: Check the number and names of files found
message("Number of genome files found: ", length(genomes_files), "\n")
if (length(genomes_files) == 0) {
  stop("No genome files found in the specified directory.")
}

for (file in genomes_files) {
  message("Found genome file: ", file, "\n")
}

message("Number of cores: ", parallel::detectCores(), "\n")

tictoc::tic.clearlog()
tictoc::tic(paste0("Running microtrait for ", length(genomes_files)))

ncores = floor(parallel::detectCores() * 0.9)

# Error handling for extract.traits.parallel
microtrait_results = tryCatch({
  extract.traits.parallel(genomes_files, dirname(genomes_files), ncores = ncores)
}, error = function(e) {
  message("Error in extract.traits.parallel: ", e$message, "\n")
  NULL
})

tictoc::toc(log = TRUE)

# Debugging: Check if microtrait results are obtained
if (is.null(microtrait_results) || length(microtrait_results) == 0) {
  stop("Microtrait extraction failed. No results obtained.")
}
print(microtrait_results)

# Error handling for unlist and parallel::mclapply
rds_files = tryCatch({
  unlist(parallel::mclapply(microtrait_results, "[[", "rds_file", mc.cores = ncores))
}, error = function(e) {
  message("Error in extracting rds files: ", e$message, "\n")
  NULL
})

# Debugging: Check if rds_files are obtained
message("Number of rds files: ", length(rds_files), "\n")
if (is.null(rds_files) || length(rds_files) == 0) {
  stop("No rds files found in the microtrait results.")
}

rds_files = list.files(file.path(genomes_dir), full.names = T, recursive = F, pattern = ".microtrait.rds$")
print(rds_files)

# Error handling for make.genomeset.results
genomeset_results = tryCatch({
  make.genomeset.results(rds_files = rds_files,
                         ids = sub(".microtrait.rds", "", basename(rds_files)),
                         ncores = ncores)
}, error = function(e) {
  message("Error in make.genomeset.results: ", e$message, "\n")
  NULL
})

if (is.null(genomeset_results)) {
  stop("Failed to create genomeset results.")
}

# saving genome set results
saveRDS(genomeset_results, file=snakemake@output[["rds"]])

message("Microtrait analysis completed successfully.\n")

