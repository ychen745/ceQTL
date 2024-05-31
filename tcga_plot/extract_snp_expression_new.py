import os
import argparse

def extract_snp_expression(infile, f_snp01, f_snp10, f_expr01, f_expr10, outfolder):
	if outfolder not in os.listdir(os.getcwd()):
		os.mkdir(os.path.join(os.getcwd(), outfolder))
	# f = open("/research/labs/pharmacology/junwenwang/data/yanxi/snp_list.txt")
	f = open(infile)
	f.readline()
	sample_list_01 = list()
	sample_list_10 = list()
	enhancer_set_01 = set()
	enhancer_set_10 = set()
	target_set_01 = set()
	target_set_10 = set()
	snp_set_01 = set()
	snp_set_10 = set()
	dict_snp_01 = dict()
	dict_snp_10 = dict()

	dict_enhancer_01 = dict()
	dict_enhancer_10 = dict()

	dict_target_01 = dict()
	dict_target_10 = dict()

	for line in f:
		linelist = line[:-1].split("\t")
		# print(len(linelist))
		rsid = linelist[0]
		enhancer = linelist[1]
		target = linelist[2]
		ref = linelist[3]
		alt = linelist[4]
		code = linelist[5]

		if code == "01":
			snp_set_01.add(rsid)
			enhancer_set_01.add(enhancer)
			target_set_01.add(target)
		elif code == "10":
			snp_set_10.add(rsid)
			enhancer_set_10.add(enhancer)
			target_set_10.add(target)
	f.close()

	# f_snp_01 = open("/research/labs/pharmacology/junwenwang/data/tcga/LUNG_luad/GENDER_01/TCGA_snp_398394_514.txt")
	# f_snp_10 = open("/research/labs/pharmacology/junwenwang/data/tcga/LUNG_luad/GENDER_10/TCGA_snp_609321_413.txt")
	# f_expr_01 = open("/research/labs/pharmacology/junwenwang/data/tcga/LUNG_luad/GENDER_01/TCGA_expr_13537_514.txt")
	# f_expr_10 = open("/research/labs/pharmacology/junwenwang/data/tcga/LUNG_luad/GENDER_10/TCGA_expr_13549_413.txt")

	if f_snp01 != "":
		f_snp_01 = open(f_snp01)
		for line in f_snp_01:
			linelist = line[:-1].split("\t")
			snp = linelist[0]
			if snp in snp_set_01:
				dict_snp_01[snp] = linelist[1:]
	if f_snp10 != "":
		f_snp_10 = open(f_snp10)
		for line in f_snp_10:
			linelist = line[:-1].split("\t")
			snp = linelist[0]
			if snp in snp_set_10:
				dict_snp_10[snp] = linelist[1:]
	if f_expr01 != "":
		f_expr_01 = open(f_expr01)
		sample_list_01 = f_expr_01.readline()[:-1].split("\t")
		for line in f_expr_01:
			linelist = line[:-1].split("\t")
			gene = linelist[0]
			if gene in enhancer_set_01:
				dict_enhancer_01[gene] = linelist[1:]
			elif gene in target_set_01:
				dict_target_01[gene] = linelist[1:]
	if f_expr10 != "":
		f_expr_10 = open(f_expr10)
		sample_list_10 = f_expr_10.readline()[:-1].split("\t")
		for line in f_expr_10:
			linelist = line[:-1].split("\t")
			gene = linelist[0]
			if gene in enhancer_set_10:
				dict_enhancer_10[gene] = linelist[1:]
			elif gene in target_set_10:
				dict_target_10[gene] = linelist[1:]

	# f = open("/research/labs/pharmacology/junwenwang/data/yanxi/snp_list.txt")
	f = open(infile)
	f.readline()
	for line in f:
		linelist = line[:-1].split("\t")
		rsid = linelist[0]
		enhancer = linelist[1]
		target = linelist[2]
		ref = linelist[3]
		alt = linelist[4]
		code = linelist[5]
		if outfolder[-1] == "/":
			outfolder = outfolder[:-1]
		fout = open(outfolder + '/' + code + "-" + target + "-" + rsid + "-" + enhancer + ".txt", "w")
		# fout = open("/research/labs/pharmacology/junwenwang/data/yanxi/plt_sample/" + code + "-" + target + "-" + rsid + "-" + enhancer + ".txt", "w")
		fout.write("sample\tgenotype\t" + enhancer + "\t" + target + "\n")
		if code == "01":
			enhancer_list = dict_enhancer_01[enhancer]
			target_list = dict_target_01[target]
			snp_list = dict_snp_01[rsid]
			for i in range(min(len(sample_list_01), len(snp_list))):
				if snp_list[i] == "0":
					genotype = ref + ref
				elif snp_list[i] == "1":
					genotype = ref + alt
				elif snp_list[i] == "2":
					genotype = alt + alt
				else:
					genotype = "NA"
				fout.write("\t".join([sample_list_01[i], genotype, enhancer_list[i], target_list[i]]) + "\n")
		elif code == "10":
			enhancer_list = dict_enhancer_10[enhancer]
			target_list = dict_target_10[target]
			snp_list = dict_snp_10[rsid]
			for i in range(min(len(sample_list_10), len(snp_list))):
				if snp_list[i] == "0":
					genotype = ref + ref
				elif snp_list[i] == "1":
					genotype = ref + alt
				elif snp_list[i] == "2":
					genotype = alt + alt
				else:
					genotype = "NA"
				fout.write("\t".join([sample_list_10[i], genotype, enhancer_list[i], target_list[i]]) + "\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--infile', dest='infile')
	parser.add_argument('--snp01', dest='snp01', default='')
	parser.add_argument('--expr01', dest='expr01', default='')
	parser.add_argument('--snp10', dest='snp10', default='')
	parser.add_argument('--expr10', dest='expr10', default='')
	parser.add_argument('--out', dest='out', default='./')

	args = parser.parse_args()
	extract_snp_expression(args.infile, args.snp01, args.snp10, args.expr01, args.expr10, args.out)