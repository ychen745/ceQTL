with open('ap.txt') as f:
	with open('ap_out.csv', 'w') as fout:
		f.readline()
		fout.write(','.join(['gene', 'tf', 'enhancer', 'snpid', 'TF-TG', 'enhancers']) + '\n')
		for line in f:
			linelist = line.split('\n')[0].split('\t')
			gene = linelist[0]
			tf = linelist[1]
			enhancer = ';'.join(linelist[2].split(','))
			snpid = linelist[3]
			tftg = linelist[4]
			enhancers_ls = []
			for ele in linelist[5:]:
				if ele != '':
					enhancers_ls.append(ele)
			enhancers = ';'.join(enhancers_ls)
			fout.write(','.join([gene, tf, enhancer, snpid, tftg, enhancers]) + '\n')