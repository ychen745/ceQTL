import pandas as pd

fsnp = 'data/hnsc/snp_mtx.txt'
output = 'MAF.txt'

# Calculate MAF

df_snp = pd.read_table(fsnp, index_col=0)
df_snp = df_snp.astype(pd.Int32Dtype(), copy=False)

with open(output, 'w') as f:
    f.write('SNP\tMAF\n')
    for i in range(len(df_snp.index)):
        row = df_snp.iloc[i, :]
        num_0 = len(row[row == 0])
        num_1 = len(row[row == 1])
        num_2 = len(row[row == 2])

        MAF = (num_1 + 2 * num_2) / (len(row) * 2)
        
        f.write(df_snp.index[i] + '\t' + str(MAF) + '\n')