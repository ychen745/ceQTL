library(ggplot2)
library(tidyverse)
library(viridis)
library(ggpubr)
library(readxl)

tcga_vio_plot <- function(filename, snp, gene, ref, alt, out, model = "default") {
  df <- read.csv(filename, sep = "\t")

  y_val <- enquo(gene)
  pdf(file.path(out, paste0(snp, "-", gene, ".pdf")))
  plt <- ggplot(df, aes(x = genotype, y = !!y_val, fill = genotype)) +
    geom_violin(trim = FALSE) +
    stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, color = "red", fill = "red") +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = paste("Violin Plot - ", gene),
         x = "Genotype",
         y = "Gene Expression") +
    theme_minimal()
  print(plt)
  dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
input_folder <- as.character(args[1])
output_folder <- as.character(args[2])
model <- as.character(args[3])

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

for (fgene in list.files(input_folder)){
  filename <- file.path(input_folder, fgene)
  fname <- unlist(strsplit(fgene, ".txt"))
  snp_gene <- unlist(strsplit(fname, "-"))
  snp_id <- snp_gene[1]
  snp <- unlist(strsplit(snp_gene[1], "_"))
  gene <- snp_gene[2]
  ref <- snp[4]
  alt <- snp[5]
  tcga_vio_plot(filename, snp_id, gene, ref, alt, output_folder, model)
}

# tcga_vio_plot("vio_20240529/chr5_1315228_rs4975615_C_T-CLPTM1L.txt", "CLPTM1L", "C", "T", "vio_result_20240529")