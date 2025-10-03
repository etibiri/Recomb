#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phangorn)
  library(data.table)
  library(grDevices)
})

opt_list <- list(
  make_option(c("--aln"), type="character", help="MSA FASTA (DNA)"),
  make_option(c("--out"), type="character", help="Output base path (no extension)"),
  make_option(c("--metadata"), type="character", default=NA, help="TSV with first column=sample_id"),
  make_option(c("--color-map"), type="character", default=NA, help="CSV output mapping taxon,color,group"),
  make_option(c("--max-tips"), type="integer", default=1000),
  make_option(c("--per-species-cap"), type="integer", default=1000),
  make_option(c("--seed"), type="integer", default=1),
  make_option(c("--pdf-width"), type="double", default=10),
  make_option(c("--pdf-height"), type="double", default=8)
)
opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$aln) || is.null(opt$out)) {
  stop("Required: --aln and --out")
}

set.seed(opt$seed)

# Ensure output dir exists
outdir <- dirname(opt$out)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pdf_out <- paste0(opt$out, ".pdf")
png_out <- paste0(opt$out, ".png")
col_csv <- if (is.na(opt$`color-map`)) file.path(outdir, "taxa_colors.csv") else opt$`color-map`

message("[INFO] NeighborNet \u2192 ", opt$out)

# --- Read alignment
X <- read.FASTA(opt$aln, type="DNA")
taxa <- names(X)
if (length(taxa) == 0L) stop("MSA empty or unreadable: ", opt$aln)

# Optional cap to avoid overplotting
if (length(taxa) > opt$`max-tips`) {
  keep <- sample(seq_along(taxa), opt$`max-tips`)
  X <- X[keep]
  taxa <- names(X)
}

# --- Distance + NeighborNet
#d <- tryCatch(dist.dna(as.DNAbin(X), model = "TN93", pairwise.deletion = TRUE),
#              error = function(e) NULL)
#if (is.null(d)) d <- dist.dna(as.DNAbin(X), model = "raw", pairwise.deletion = TRUE)
#net <- neighborNet(d)
d <- tryCatch(dist.dna(X, model="TN93", pairwise.deletion=TRUE), error=function(e) NULL)
if (is.null(d)) dist_ok <- dist.dna(X, model="raw", pairwise.deletion=TRUE)
net <- neighborNet(d)

# --- Colors from metadata (optional)
tip_cols <- rep("#444444", length(taxa))
colmap <- c(ALL = "#444444")
grp_vec <- rep("ALL", length(taxa))

if (!is.na(opt$metadata) && file.exists(opt$metadata)) {
  meta <- tryCatch(fread(opt$metadata, sep = "\t", header = TRUE), error = function(e) NULL)
  if (!is.null(meta) && ncol(meta) >= 2) {
    # candidate columns in order of preference
    cand <- c("species", "host", "country", "year")
    col_hit <- cand[cand %in% names(meta)][1]
    # assume first column is sample_id
    m <- meta[match(taxa, meta[[1]]), ]
    if (!is.na(col_hit)) {
      grp_vec <- as.character(m[[col_hit]])
      grp_vec[is.na(grp_vec) | grp_vec == ""] <- "NA"
      ug <- unique(grp_vec)
      pal <- hcl.colors(length(ug), "Dark3", rev = FALSE)
      colmap <- setNames(pal, ug)
      tip_cols <- unname(colmap[grp_vec])
    }
  }
}

# --- Write color map
col_dt <- data.table(taxon = taxa, color = tip_cols, group = grp_vec)
fwrite(col_dt, file = col_csv, sep = ",", quote = FALSE)

# --- Draw PDF
pdf(pdf_out, width = opt$`pdf-width`, height = opt$`pdf-height`, family = "Helvetica")
plot(net, "2D", show.tip.label = FALSE, edge.color = "#999999", tip.color = tip_cols)
legend("topleft", legend = names(colmap), col = unname(colmap), pch = 16, bty = "n", cex = 0.8)
dev.off()

# --- Draw PNG (96 dpi)
png(png_out, width = as.integer(96 * opt$`pdf-width`), height = as.integer(96 * opt$`pdf-height`), res = 96)
plot(net, "2D", show.tip.label = FALSE, edge.color = "#999999", tip.color = tip_cols)
legend("topleft", legend = names(colmap), col = unname(colmap), pch = 16, bty = "n", cex = 0.8)
dev.off()

message("[OK] Wrote: ", pdf_out, " | ", png_out, " | ", col_csv)
