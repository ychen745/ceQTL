import pandas as pd
import scipy
from statsmodels.stats.multitest import multipletests

def anova_all(fsnp, fexpr, fpair, output):
    df_snp = pd.read_table(fsnp, index_col=0)
    df_expr = pd.read_table(fexpr, index_col=0)

    df_snp = df_snp.astype(pd.Int32Dtype(), copy=False)
    df_expr = df_expr.astype(pd.Float64Dtype(), copy=False)

    intersect = list(set(df_snp.columns).intersection(df_expr.columns))

    df_snp = df_snp[intersect]
    df_expr = df_expr[intersect]

    with open(fpair) as f:
        with open(output, 'w') as fout:
            f.readline()
            fout.write(','.join(['gene', 'tf', 'snp_id', 'TF-TG?', 'p_gene', 'p_tf']) + '\n')
            for line in f:
                linelist = line.split('\n')[0].split(',')
                gene = linelist[0]
                tf = linelist[1]
                enhancer = linelist[2]
                snpid = linelist[3]
                tftg = linelist[4]
                if snpid not in df_snp.index:
                    continue
                if gene not in df_expr.index or tf not in df_expr.index:
                    continue
                genotypes = df_snp.loc[snpid]
                group0 = genotypes[genotypes == 0].index
                group1 = genotypes[genotypes == 1].index
                group2 = genotypes[genotypes == 2].index

                gene0 = df_expr.loc[gene, group0]
                gene1 = df_expr.loc[gene, group1]
                gene2 = df_expr.loc[gene, group2]

                tf0 = df_expr.loc[tf, group0]
                tf1 = df_expr.loc[tf, group1]
                tf2 = df_expr.loc[tf, group2]

                stat_gene, p_gene = scipy.stats.f_oneway(gene0, gene1, gene2)
                stat_tf, p_tf = scipy.stats.f_oneway(tf0, tf1, tf2)

                fout.write(','.join([gene, tf, snpid, tftg, str(p_gene), str(p_tf)]) + '\n')

def p_adj(fres, output):
    p_gene = []
    p_tf = []

    with open(fres) as f:
        f.readline()
        for line in f:
            linelist = line.split('\n')[0].split(',')
            p_gene.append(float(linelist[4]))
            p_tf.append(float(linelist[5]))

    _, p_bonf_gene, _, _ = multipletests(p_gene, alpha=0.05, method='bonferroni')
    _, p_bh_gene, _, _ = multipletests(p_gene, alpha=0.05, method='fdr_bh')

    _, p_bonf_tf, _, _ = multipletests(p_tf, alpha=0.05, method='bonferroni')
    _, p_bh_tf, _, _ = multipletests(p_tf, alpha=0.05, method='fdr_bh')

    idx = 0

    with open(output, 'w') as fout:
        fout.write(','.join(['gene', 'tf', 'snp_id', 'TF-TG?', 'p_gene', 'p_gene_bonf', 'p_gene_bh', 'p_tf', 'p_tf_bonf', 'p_tf_bh']) + '\n')
        with open(fres) as f:
            f.readline()
            for line in f:
                linelist = line.split('\n')[0].split(',')
                gene = linelist[0]
                tf = linelist[1]
                snp_id = linelist[2]
                tftg = linelist[3]
                p_gene = linelist[4]
                p_tf = linelist[5]
                fout.write(','.join([gene, tf, snp_id, tftg, p_gene, str(p_bonf_gene[idx]), str(p_bh_gene[idx]), p_tf, str(p_bonf_tf[idx]), str(p_bh_tf[idx])]) + '\n')
                idx += 1

if __name__ == '__main__':
    fsnp = 'data/hnsc/snp_mtx.txt'
    fexpr = 'data/hnsc/expr_mtx.txt'
    fpair = 'ap.csv'
    output = 'result.csv'

    fres = 'hnsc_anova_20240529.csv'
    fout = 'hnsc_anova_20240529_adj.csv'

    # anova_all(fsnp, fexpr, fpair, output)
    # p_adj(fres, fout)



