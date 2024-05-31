library("tidyverse")
library("MatrixEQTL")

base.dir = find.package("MatrixEQTL")

fsnp = file.path("data", "hnsc", "snp_mtx.txt")
fgene = file.path("data", "hnsc", "expr_mtx.txt")
fsnp_loc = file.path("data", "hnsc", "snpsloc.txt")
fgene_loc = file.path("data", "reference", "genepos_hg19.txt")

snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters = "NA";
snps$fileSkipRows = 1;
snps$fileSkipColumns = 1;
snps$fileSliceSize = 2000;
snps$LoadFile(fsnp);

gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA";
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 2000;
gene$LoadFile(fgene);

snpspos = read.table(fsnp_loc, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(fgene_loc, header = TRUE, stringsAsFactors = FALSE);

fout_cis = file.path("results", "hnsc", "cis.txt");
fout_trans = file.path("results", "hnsc", "trans.txt");

p_threshold_cis = 1e-4;
p_threshold_trans = 1e-4;
cisDist = 1e6;

useModel = modelLINEAR;

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  output_file_name = NULL,
  pvOutputThreshold = 0,
  useModel = useModel,
  verbose = TRUE,
  output_file_name.cis = fout_cis,
  pvOutputThreshold.cis = p_threshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

show(me$cis$eqtls)