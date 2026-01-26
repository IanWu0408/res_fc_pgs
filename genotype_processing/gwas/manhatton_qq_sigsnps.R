
#!/usr/bin/env Rscript
# ============================================================
# GWAS visualization & SNP-to-gene mapping (single-file)
#
# This script:
#   1) Reads GWAS summary statistics
#   2) Generates Manhattan and QQ plots
#   3) Computes genomic inflation factor (lambda)
#   4) Maps significant SNPs to nearby genes (±window)
#   5) Exports SNP–gene tables
#   6) Generates volcano plots at multiple p-value thresholds
#
# Author: Tung-Yen Wu
# ============================================================


# ============================================================
# 0) USER CONFIG — EDIT THIS SECTION ONLY
# ============================================================

CFG <- list(
  
  ## ---- Input ----
  sumstats_path = "/PATH/TO/gwas.sumstats.gz",
  
  ## ---- Output directory ----
  out_dir = "/PATH/TO/output_dir",
  
  ## ---- Plot parameters ----
  genomewide_p = 5e-8,
  suggestive_p = 1e-5,
  volcano_beta_xlim = c(-0.1, 0.1),
  
  ## ---- SNP → gene mapping ----
  gene_window_bp = 5000,     # ± window for gene assignment
  p_thresh_gene  = 1e-5      # SNPs used for gene mapping
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1) Packages
# ============================================================

suppressPackageStartupMessages({
  pkgs <- c(
    "data.table", "dplyr", "ggplot2",
    "qqman", "RColorBrewer", "R.utils",
    "BiocManager", "GenomicRanges", "biomaRt"
  )
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (p %in% c("GenomicRanges", "biomaRt")) {
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p)
      }
    }
  }
  
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(qqman)
  library(RColorBrewer)
  library(R.utils)
  library(GenomicRanges)
  library(biomaRt)
})


# ============================================================
# 2) Load & clean GWAS summary statistics
# ============================================================

dat <- fread(CFG$sumstats_path)
setDT(dat)

# Handle chromosome X
dat[CHR == "X", CHR := "23"]

manhattan_data <- dat[, .(
  SNP  = SNP,
  CHR  = as.numeric(CHR),
  BP   = as.numeric(BP),
  P    = as.numeric(P),
  BETA = as.numeric(BETA),
  F    = FREQ,
  SE   = SE
)]

manhattan_data <- manhattan_data[!is.na(P) & P > 0]
manhattan_data[, logp := -log10(P)]


# ============================================================
# 3) Manhattan plot
# ============================================================

chr_len <- manhattan_data[, .(max_bp = max(BP)), by = CHR]
chr_len[, offset := cumsum(max_bp) - max_bp]
manhattan_data <- merge(manhattan_data, chr_len, by = "CHR")
manhattan_data[, cum_pos := BP + offset]

axis_df <- manhattan_data[, .(
  center = (min(cum_pos) + max(cum_pos)) / 2
), by = CHR]

colors <- brewer.pal(8, "Dark2")

p_manhattan <- ggplot(manhattan_data,
                      aes(x = cum_pos, y = logp, color = factor(CHR))) +
  geom_point(size = 0.6) +
  geom_hline(yintercept = -log10(CFG$suggestive_p),
             linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(CFG$genomewide_p),
             linetype = "dashed", color = "black") +
  scale_color_manual(values = rep(colors, length.out = length(unique(manhattan_data$CHR)))) +
  scale_x_continuous(breaks = axis_df$center,
                     labels = axis_df$CHR) +
  labs(x = "Chromosome", y = "-log10(p-value)",
       title = "GWAS Manhattan Plot") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(file.path(CFG$out_dir, "manhattan.png"),
       p_manhattan, width = 12, height = 6)


# ============================================================
# 4) QQ plot & lambda
# ============================================================

png(file.path(CFG$out_dir, "qq_plot.png"),
    width = 600, height = 600, res = 120)
qq(manhattan_data$P, main = "GWAS QQ Plot")
dev.off()

chisq  <- qchisq(1 - manhattan_data$P, 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

writeLines(paste("Genomic inflation factor (lambda):", round(lambda, 4)),
           file.path(CFG$out_dir, "lambda.txt"))


# ============================================================
# 5) SNP → gene mapping (±window)
# ============================================================

sig_snps <- manhattan_data[P <= CFG$p_thresh_gene]

snps_gr <- GRanges(
  seqnames = sig_snps$CHR,
  ranges   = IRanges(
    start = sig_snps$BP - CFG$gene_window_bp,
    end   = sig_snps$BP + CFG$gene_window_bp
  ),
  SNP = sig_snps$SNP
)

mart <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "useast"
)

genes_df <- getBM(
  attributes = c("chromosome_name","start_position",
                 "end_position","hgnc_symbol"),
  filters    = "chromosome_name",
  values     = as.character(unique(sig_snps$CHR)),
  mart       = mart
)

genes_df <- genes_df[hgnc_symbol != ""]
genes_gr <- GRanges(
  seqnames = genes_df$chromosome_name,
  ranges   = IRanges(
    start = genes_df$start_position,
    end   = genes_df$end_position
  ),
  gene = genes_df$hgnc_symbol
)

hits <- findOverlaps(snps_gr, genes_gr)

snp2gene <- data.table(
  SNP  = sig_snps$SNP[queryHits(hits)],
  Gene = genes_gr$gene[subjectHits(hits)],
  P    = sig_snps$P[queryHits(hits)],
  BETA = sig_snps$BETA[queryHits(hits)],
  F    = sig_snps$F[queryHits(hits)],
  SE   = sig_snps$SE[queryHits(hits)]
)

fwrite(unique(snp2gene),
       file.path(CFG$out_dir, "snp2gene_p1e5.csv"))


# ============================================================
# 6) Volcano plots
# ============================================================

p_thresholds <- c(
  "5e-8" = CFG$genomewide_p,
  "1e-5" = CFG$suggestive_p
)

for (nm in names(p_thresholds)) {
  
  p_thr <- p_thresholds[nm]
  
  df <- manhattan_data %>%
    mutate(sig = ifelse(P < p_thr, "sig", "not sig"))
  
  p_volcano <- ggplot(df, aes(x = BETA, y = logp, color = sig)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("not sig" = "gray70", "sig" = "red")) +
    geom_hline(yintercept = -log10(p_thr),
               linetype = "dashed", color = "blue") +
    coord_cartesian(xlim = CFG$volcano_beta_xlim) +
    labs(
      x = expression(beta),
      y = expression(-log[10](p)),
      title = paste0("Volcano plot (p < ", nm, ")")
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  ggsave(
    file.path(CFG$out_dir, paste0("volcano_p", nm, ".png")),
    p_volcano,
    width = 12,
    height = 6
  )
}

cat("[done] GWAS visualization & SNP-to-gene mapping completed\n")