#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ape); library(phangorn) })
args <- commandArgs(trailingOnly=TRUE)

getArg <- function(flag, default=NA_character_) {
  i <- match(flag, args); if (is.na(i) || i>=length(args)) return(default); args[i+1]
}
aln <- getArg("--aln"); out <- getArg("--out")
max_tips <- as.integer(getArg("--max-tips", "500"))
seed <- as.integer(getArg("--seed", "42"))

if (is.na(aln) || is.na(out) || nchar(aln)==0 || nchar(out)==0) {
  cat("[FATAL] Missing --aln/--out\n", file=stderr()); quit(status=2)
}
odir <- dirname(out); if (!dir.exists(odir)) dir.create(odir, recursive=TRUE, showWarnings=FALSE)

safe_plot <- function(expr, title_fallback) {
  pdf(out, width=8, height=6)
  ok <- TRUE
  tryCatch({ eval.parent(substitute(expr)) }, error=function(e){
    ok <<- FALSE; plot.new(); title(main=paste(title_fallback, "\n", conditionMessage(e)))
  })
  dev.off(); ok
}

ok <- safe_plot({
  X <- read.dna(aln, format="fasta")
  n <- if (is.null(X)) 0L else nrow(X)
  if (n < 2) stop("Alignment has <2 sequences")
  if (!is.na(max_tips) && n > max_tips) {
    set.seed(seed)
    keep <- sample(seq_len(n), size=max_tips, replace=FALSE)
    X <- X[keep, ]
  }
  d <- dist.dna(X, model="raw", pairwise.deletion=TRUE)
  net <- neighborNet(d)
  plot(net, "2D", show.tip.label=FALSE)
}, sprintf("NeighborNet placeholder for %s", basename(aln)))

cat(if (ok) sprintf("[INFO] NeighborNet → %s\n", out)
    else sprintf("[WARN] NeighborNet placeholder → %s\n", out))
