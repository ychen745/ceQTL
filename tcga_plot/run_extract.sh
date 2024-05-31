INPUT='plot_20240529.tsv'
SNP='../data/hnsc/snp_mtx.txt'
EXPR='../data/hnsc/expr_mtx.txt'
OUTPUT='plt_20240529'

python extract_snp_expression_new.py --infile $INPUT --snp01 $SNP --snp10 $SNP --expr01 $EXPR --expr10 $EXPR --out $OUTPUT
