INPUT='vio_20240529.csv'
SNP='../data/hnsc/snp_mtx.txt'
EXPR='../data/hnsc/expr_mtx.txt'
OUTPUT='vio_20240529'

python extract_snp_expr_only.py --infile ${INPUT} --snp ${SNP} --expr ${EXPR} --out ${OUTPUT}
