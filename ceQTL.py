import numpy as np
import pandas as pd

def ceQTL(fsnp, fexpr, fgene_info, ftf_info, fcov=None, ext_dist=10000):
    df_snp = pd.read_table(fsnp, index_col=0)
    df_expr = pd.read_table(fexpr)
    # df_ftf = pd.read_table(ftf_info)

    df_snp = df_snp.astype(pd.Int32Dtype(), copy=False)
    df_expr = df_expr.astype(pd.Float64Dtype(), copy=False)

    index_list = [ele.split('_') for ele in df_snp.index]
    ref_list = [ele[-2] for ele in index_list]
    alt_list = [ele[-1] for ele in index_list]

    df_snp['ref'] = ref_list
    df_snp['alt'] = alt_list

    print(df_snp.columns)


if __name__ == '__main__':
    fsnp = 'data/snp_mtx.txt'
    fexpr = 'data/expr_mtx.txt'
    ceQTL(fsnp, fexpr, None, None)