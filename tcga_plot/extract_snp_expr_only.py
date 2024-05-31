import os
import pandas as pd
import argparse

def extract_snp_expression(infile, fsnp, fexpr, outfolder):
	if outfolder not in os.listdir(os.getcwd()):
		os.mkdir(os.path.join(os.getcwd(), outfolder))

	df_snp = pd.read_table(fsnp, index_col=0)
	df_expr = pd.read_table(fexpr, index_col=0)

	df_snp = df_snp.astype(pd.Int32Dtype(), copy=False)
	df_expr = df_expr.astype(pd.Float64Dtype(), copy=False)

	intersect = list(set(df_snp.columns).intersection(df_expr.columns))

	df_snp = df_snp[intersect]
	df_expr = df_expr[intersect]


	with open(infile) as f:
		f.readline()
		for line in f:
			linelist = line.split('\n')[0].split(',')
			snp = linelist[0]
			gene = linelist[1]
			ref = snp.split('_')[-2]
			alt = snp.split('_')[-1]

			geno_dict = {0: ref + ref, 1: ref + alt, 2: alt + alt}

			snp_sub = df_snp.loc[snp, :]
			expr_sub = df_expr.loc[gene, :]

			with open(os.path.join(outfolder, snp + '-' + gene + '.txt'), 'w') as fout:
				fout.write('\t'.join(['sample', 'genotype', gene]) + '\n')
				for i in range(len(snp_sub)):
					fout.write('\t'.join([snp_sub.index[i], geno_dict[snp_sub.iloc[i]], str(expr_sub.loc[snp_sub.index[i]])]) + '\n')

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--infile')
	parser.add_argument('--snp')
	parser.add_argument('--expr')
	parser.add_argument('--out')

	args = parser.parse_args()
	extract_snp_expression(args.infile, args.snp, args.expr, args.out)