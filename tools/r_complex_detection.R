#!/usr/bin/env Rscript
# Run gapseq's complex_detection() on a JSON-encoded input batch and emit a
# JSON result, so the Rust integration test can diff against the canonical
# R implementation.
#
# Input (stdin): JSON array of objects [{ "rxn": "...", "descs": ["...","..."] }, ...]
# Output (stdout): JSON array of arrays matching the input order.

suppressMessages(library(data.table))
suppressMessages(library(stringr))

gapseq_root <- Sys.getenv("GAPSEQ_ROOT", "/home/leechuck/Public/software/gapseq")
subunit_dict_path <- file.path(gapseq_root, "dat/complex_subunit_dict.tsv")

# Load the real function from gapseq's source.
# complex_detection.R has a top-level `source` statement for the dict at its
# tail — but running it directly also defines `complex_detection()` at the top
# level. We source it with the script.dir pointing to the src/ folder so the
# file's own `script.dir` resolution works.
script.dir <- file.path(gapseq_root, "src")
source(file.path(script.dir, "complex_detection.R"))

# After sourcing, the function name `complex_detection` is defined.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
    inpath <- args[1]
    outpath <- if (length(args) >= 2) args[2] else "-"
    raw <- paste(readLines(inpath, warn = FALSE), collapse = "\n")
} else {
    raw <- paste(readLines("stdin", warn = FALSE), collapse = "\n")
    outpath <- "-"
}

# Minimal JSON parser substitute using eval() — our input is always simple.
# Avoid pulling `jsonlite` in since we don't want the dependency.
# Instead, use a plain text format: each line is `rxn<TAB>desc1|||desc2|||...`.
lines <- strsplit(raw, "\n", fixed = TRUE)[[1]]
lines <- lines[nchar(lines) > 0]

results <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    rxn <- parts[1]
    descs <- if (length(parts) >= 2 && nchar(parts[2]) > 0) {
        strsplit(parts[2], "|||", fixed = TRUE)[[1]]
    } else character(0)
    res <- complex_detection(descs, rxn)
    res[is.na(res)] <- ""
    paste(res, collapse = "|||")
})

out <- paste(results, collapse = "\n")
if (outpath == "-") {
    cat(out, "\n", sep = "")
} else {
    writeLines(out, outpath)
}
