# input file tsv 6 columns:
rsid, enhancer, target, ref, alt, code

# extract expression
python --infile snp_list.txt --snp01 snp_mtx_01 --expr01 expr_mtx_01 --snp10 snp_mtx_10 --expr10 expr_mtx_10 --out output_folder

# draw plots, creates a subfolder 'plots'
Rscript tcga_cor_plot.r snp_list.txt
