DEBUG = TRUE

##################################################################   find snp 

ceQTL = function(genename, fgeno, fexpr, fgene_loc, ftf_info, fcov = '', extend_dis = 10000){
  if(!require("tidyverse")){
    install.packages("tidyverse")
  }
  if(!"tidyverse" %in% (.packages())){
    library(tidyverse)
  }
  
  ### debugging ###
  
  # if(DEBUG){
  #   fgeno = "/research/labs/pharmacology/junwenwang/data/yanxi/tcga_test/snp_mtx.txt"
  #   fexpr = "/research/labs/pharmacology/junwenwang/data/yanxi/tcga_test/expr_mtx.txt"
  #   fgene_loc = "/research/labs/pharmacology/junwenwang/data/yanxi/reference/TCGA_gene_loc_strand_hg19.txt"
  #   ftf_info = "/research/labs/pharmacology/junwenwang/data/genecard/all/info/AACS.txt"
  # }
  
  ### end debugging ###
  
  snp = read.csv(fgeno, header = T, sep="\t", row.names = 1)
  snp = checkSnpName(snp)
  expr = read.csv(fexpr, header = T, sep="\t", row.names = 1)
  genotype = extractGenotype(snp)
  genotype = t(data.frame(as.list(genotype)))
  colnames(genotype) = c("genotype")
  genotype = as.data.frame(cbind(snpid = rownames(genotype), genotype))
  gene_loc = read.table(fgene_loc, header = T, sep = '\t')
  snp_loc = extractSnpLoc(snp)
  snp_loc = t(data.frame(as.list(snp_loc)))
  colnames(snp_loc) = c("snploc")
  snp_loc = as.data.frame(cbind(rsid = rownames(snp_loc), snp_loc))
  snp_loc = snp_loc %>%
    separate(snploc, into = c("chr", "position"), sep = ":") %>%
    select(c("rsid", "chr", "position"))
  
  gene_list = read.table(ftf_info, header = T, sep = '\t')
  gene_list = select(gene_list, c("gh_id", "location_hg38", "tss_distance", "tfbs"))
  # gene_list = gene_list[,c(1,5,4,7)]
  colnames(gene_list) = c("Enhancer_id","Loc","tss_dis","TFs" )
  gene_list$TFs = as.character(gene_list$TFs)
  
  result = findsnp(genename,gene_list,gene_loc,snp_loc,extend_dis)   # find snp
  
  if(nrow(result) < 1) {
    print("Error: No overlapping SNP found.")
    stop()
  }
  
  snp_gene = case_inputfile(expr,snp,result,gene_list,genename)
  
  if(nrow(snp_gene) < 1) {
    print("Error: failed to write case file.")
    stop()
  }
  
  case = read.table(paste0(genename,'_case_snp_gene.xls'),header = F, sep = '\t')
  
  result_anova = CalAnovaPccs(snp_gene,expr,snp)   # caculat anova and pccs
  snp_inputfile(snp, snp_gene, case, genename)  #caculate snp file for chowtest
  expr_inputfile(expr,genename)   # create expr file for chowtest
  
  ### caculate chowtest ###
  
  infile   =  paste0(genename,'_expr_all.txt')
  gfile    =  paste0(genename,'_snp_all.txt')
  cfile    =  paste0(genename,'_case_snp_gene.xls')
  outfile  =  paste0(genename,'_chow_result_ALL.txt')
  covfile = fcov
  
  pan_chowtest (infile,gfile,cfile,outfile,covfile,pcut)  # chowtes all
  
  ### end chowtest ###
  
  chow_group= read.table(outfile, header = T, sep = '\t')
  
  final_all = caculateALL(chow_group,result)
  final_all = merge(final_all,result_anova,by=c('enhancer','snpid','tfgene'))
  final_all = merge(final_all,genotype,by='snpid')
  final_all = final_all[order(abs(final_all$Pval),decreasing = F),]
  final_all = final_all[,c(4,2,1,3,24,5,15:23,6:14)]
  colnames(final_all)[6] = 'chow_pval'
  df = which(duplicated(paste0(final_all$snpid,final_all$tfgene)))
  if(length(df)>0) {
    final_dup = final_all[df,2:4]
    final_dup$uid = paste0(final_dup$snpid,final_dup$tfgene)
    final_all = final_all[-df,]
    for( i in 1:length(df))
    {
      tmp = which(paste0(final_all$snpid,final_all$tfgene)== final_dup[i,]$uid)
      final_all[tmp,]$enhancer=paste0(final_all[tmp,]$enhancer,',',final_dup[i,]$enhancer)
    }
  }
  
  outpath  = './results/'
  
  filename  =  paste0(outpath, genename, '_ALL.txt')
  
  if(nrow(final_all) > 0){
    write.table(final_all,filename,row.names = F,sep = '\t')
  }
}

checkSnpName = function(snp_mtx){
  str = rownames(snp_mtx)[1]
  chrstr = unlist(strsplit(str, "_"))
  if(length(chrstr) != 5){
    chrstr = unlist(strsplit(str, "\\."))
    if(length(chrstr) != 5){
      print(str)
      print("Error: unrecognized SNP ID format.")
      stop()
    }
  } else{
    rownames(snp_mtx) = sapply(rownames(snp_mtx), fixSnpName)
  }
  return(snp_mtx)
}

fixSnpName = function(str){
  return(paste0(unlist(strsplit(str, "_")), collapse = "."))
}

extractGenotype = function(snp_mtx){
  genotypes = sapply(rownames(snp_mtx), genotypeFromString)
}

extractSnpLoc = function(snp){
  snp_loc = sapply(rownames(snp), snpLocFromString)
}

genotypeFromString = function(str) 
{
  chrstr = unlist(strsplit(str, "\\."))
  ref = chrstr[4]
  alt = chrstr[5]
  genotype = paste0(paste0(ref, ref), "/", paste0(ref, alt), "/", paste0(alt, alt))
  return(genotype)
}

snpLocFromString = function(str){
  chrstr = unlist(strsplit(str, "\\."))
  chm = chrstr[1]
  pos = chrstr[2]
  snp_loc = paste0(chm, ":", pos)
  return(snp_loc)
}

findsnp = function(genename,gene_list,gene_loc,snp_loc,extend_dis)   # find snp
{
  
  result<-matrix(nrow=T,ncol=12)
  result<-as.data.frame(result)
  result=na.omit(result)
  
  for( i in 1:nrow(gene_list))   # for each enhancer find snp
  {
    tmp_result<-matrix(nrow=T,ncol=10)
    tmp_result<-as.data.frame(tmp_result)
    tmp_result=na.omit(tmp_result)
    enhancer = paste(gene_list[i,]$Enhancer_id )
    snprange = paste(gene_list[i,]$Loc) 
    tss_dis  = gene_list[i,]$tss_dis 
    
    chrstr = unlist(strsplit(snprange,':'))
    chrstr.1 = unlist(strsplit(chrstr[2],'-'))
    chrstr = unlist(c(chrstr[1],chrstr.1[1],chrstr.1[2]))
    
    str.1 = as.numeric(chrstr[2])
    str.2 = as.numeric(chrstr[3])
    
    chrstr[2] =  as.integer(chrstr[2])-extend_dis
    chrstr[3] = as.integer(chrstr[3])+extend_dis
    #new loc
    genesnp = snp_loc[which(snp_loc$chr ==chrstr[1]),]
    
    # print(nrow(genesnp))
    
    genesnp$position = as.integer(genesnp$position)
    genesnp = genesnp[which(genesnp$position >= as.integer(chrstr[2]) & genesnp$position <= as.integer(chrstr[3])),]   # find snp pos
    
    # print(nrow(genesnp))
    
    colnames(genesnp)[2] ='chr' 
    genename_lo = gene_loc[which(gene_loc$geneid ==genename ),]   # genename loacation

    # print(genename_lo[1,])
    
    genename_snp = merge(genename_lo,genesnp,by='chr')
    
    # print(nrow(genename_snp))
    
    if(nrow(genename_snp)>0) {
      genename_snp$tss_snp_dist<-ifelse(genename_snp$strand > 0 ,genename_snp$position - genename_snp$s1,genename_snp$s2-genename_snp$position )
      genename_snp$tss_snp_dist<-ifelse(genename_snp$s1<= genename_snp$position & genename_snp$position <= genename_snp$s2,0,genename_snp$tss_snp_dist)
      colnames(tmp_result)<-colnames(genename_snp)
      tmp_result<-rbind(tmp_result,genename_snp)
      tmp_result$enhancer = enhancer
      tmp_result$enhancer_loc = snprange
      tmp_result$tss_dis = tss_dis
      strand = tmp_result[1,]$strand
      
      if(strand > 0) 
        tmp_result$enhancer_gene_dis<-ifelse(as.integer(str.2)< tmp_result$s1 , as.integer(str.2) -tmp_result$s1, as.integer(str.1) - tmp_result$s2)
      if(strand < 0)  
        tmp_result$enhancer_gene_dis<-ifelse(as.integer(str.2)< tmp_result$s1 , -(as.integer(str.2) -tmp_result$s1), -(as.integer(str.1) - tmp_result$s2))
      tmp_result$enhancer_gene_dis<-ifelse((tmp_result$s1 >= as.integer(str.1) & tmp_result$s1 <=as.integer(str.2)) | (tmp_result$s2 >= as.integer(str.1) & tmp_result$s2 <=as.integer(str.2)),0,tmp_result$enhancer_gene_dis )
      
      if(strand > 0) 
        tmp_result$enhancer_snp_dis = ifelse(tmp_result$position <as.integer(str.1),tmp_result$position - as.integer(str.1),tmp_result$position - as.integer(str.2) )
      if(strand <0)
        tmp_result$enhancer_snp_dis = ifelse(tmp_result$position <as.integer(str.1),-(tmp_result$position - as.integer(str.1)),-(tmp_result$position - as.integer(str.2)) )
      tmp_result$enhancer_snp_dis = ifelse(tmp_result$position  >= as.integer(str.1) & tmp_result$position <= as.integer(str.2),0,tmp_result$enhancer_snp_dis)
      
      tmp_result = tmp_result[,c(2,9:10,6,3:5,7:8,11:13)]
    }
    if(nrow(tmp_result)>0)
    { 
      colnames(result) = colnames(tmp_result)
      result=rbind(result,tmp_result)
    }
  }
  colnames(result)[5:6]=c('gene_start','gene_end')
  colnames(result)[8]=c('snp_pos')
  return(result)
}

expr_inputfile = function(expr,genename) {
  expr = t(expr)
  write.table(expr, paste0(genename, "_expr_all.txt"), row.names = T, sep = '\t', quote = F)
}

case_inputfile = function(expr,snp,result,gene_list,genename) {
  
  tf_genelist = matrix(nrow=T,ncol=2)
  tf_genelist = as.data.frame(tf_genelist)
  tf_genelist = na.omit(tf_genelist)
  
  for(i in 1:nrow(gene_list)) {     # split tfgene for each enhancer
    enhancer = gene_list[i,]$Enhancer_id
    chrstr= paste0(gene_list[i,]$TFs)
    chrstr=strsplit(chrstr, ';')
    str = unlist(chrstr)
    str = unique(str)
    tmp_genelist = matrix(nrow=length(str), ncol=2)
    tmp_genelist = as.data.frame(tmp_genelist)
    tmp_genelist[,1] = enhancer
    tmp_genelist[,2] = str
    tf_genelist = rbind(tf_genelist,tmp_genelist)
  }
  
  result.2 = result[,c(2,4)]
  colnames(tf_genelist) = c('enhancer','tfgene')
  snp_gene = merge(result.2,tf_genelist,by='enhancer')   # merge rsid and tf_gene
  snp_gene = snp_gene[which(snp_gene$tfgene %in% rownames(expr)),]  #check each snpid and geneid whether in matrix
  snp_gene = snp_gene[which(snp_gene$rsid %in% rownames(snp)),]
  df = which(snp_gene$tfgene == genename)
  if(length(df)>0)
    snp_gene = snp_gene[-df,]
  
  tmp_case = matrix(nrow=T,ncol=2)
  tmp_case =na.omit(tmp_case)
  
  for(i in 1:nrow(snp_gene))
  {
    tmp_case.1 = paste0(snp_gene[i,1],'_',snp_gene[i,2],'_',genename,'_',snp_gene[i,3])
    tmp_case.2 = paste(genename,'~',snp_gene[i,3])
    tmp_case = rbind(tmp_case,c(tmp_case.1,tmp_case.2))
  }    
  
  write.table(tmp_case,paste0(genename,'_case_snp_gene.xls'),row.names = F,col.names = F, sep = '\t',quote = F)  # for chowtest
  
  return(snp_gene)
  
}

snp_inputfile = function(snp,snp_gene,tmp_case,genename) {
  snp = snp[snp_gene[,2],]
  snp_matrix.all = t(snp)
  colnames(snp_matrix.all) = tmp_case[,1]
  write.table(snp_matrix.all,file=paste0(genename,"_snp_all.txt"),row.names = T,sep="\t",quote=F)
}

caculateALL = function(chow_group,result,p_cut=0.05) {
  # Case and p-value
  chow_group= chow_group[,c(1,7)]
  str = unlist(strsplit(as.character(chow_group$Case),'_'))
  chrstr <- matrix(str, ncol = 4, byrow = TRUE)
  chow_group$enhancer =  chrstr[,1]
  chow_group$snpid = chrstr[,2]
  chow_group$gene = chrstr[,3]
  chow_group$tfgene = chrstr[,4]
  chow_group = chow_group[,c(5,3:4,6,2)]
  
  result = result[,2:12]
  colnames(result)[3]='snpid'
  colnames(chow_group)[2]='enhancer'
  
  final = merge(chow_group,result,by=c('enhancer','snpid'))  #group chow_result and snp-result
  final = final[order(final$Pval,decreasing = F),]
  final = final[which(final$Pval < p_cut),]
  
  return(final)
}

CalAnovaPccs = function(snp_gene,expr,snp) {
  result = matrix(nrow = T,ncol = 12)
  result = na.omit(result)
  
  for( i in 1:nrow(snp_gene))
  {
    p0 = p1 = p2 = NULL
    r0 = r1 = r2 = NULL
    mean_gene_expr_group0 = mean_gene_expr_group1 =mean_gene_expr_group2 = -1
    mean_tf_expr_group0 =  mean_tf_expr_group1 =  mean_tf_expr_group2 = -1
    sd_gene_expr_group0 = sd_gene_expr_group1 = sd_gene_expr_group2 = -1
    sd_tf_expr_group0 = sd_tf_expr_group1 = sd_tf_expr_group2 = -1
    
    snp_value = snp[which(rownames(snp) %in%snp_gene[i,]$rsid),]
    
    df<-which( snp_value ==0)
    group_0<-colnames(snp_value)[df]
    df<-which(snp_value ==1 )
    group_1<-colnames(snp_value)[df]
    df<-which(snp_value ==2 )
    group_2<-colnames(snp_value)[df]
    
    gene_expr = expr[which(rownames(expr) %in% genename),]
    gene_expr_group0 = as.numeric(gene_expr[,which(colnames(gene_expr) %in% group_0)])
    gene_expr_group1 = as.numeric(gene_expr[,which(colnames(gene_expr) %in% group_1)])
    gene_expr_group2 = as.numeric(gene_expr[,which(colnames(gene_expr) %in% group_2)])
    
    tf_expr = expr[which(rownames(expr) %in% snp_gene[i,]$tfgene ),]
    
    tf_expr_group0 = as.numeric(tf_expr[,which(colnames(tf_expr) %in% group_0)] )
    tf_expr_group1 = as.numeric( tf_expr[,which(colnames(tf_expr) %in% group_1)] )
    tf_expr_group2 = as.numeric( tf_expr[,which(colnames(tf_expr) %in% group_2)] )
    
    if(length(group_0) < 2 || length(group_1) < 2 || length(group_2) < 2){
      print("group length less than 2!")
      stop()
    }
    
    # group 0
    mean_gene_expr_group0 = round(mean(gene_expr_group0),4)
    mean_tf_expr_group0 = round(mean(tf_expr_group0),4)
    sd_tf_expr_group0 = round(sd(tf_expr_group0)*sqrt((length(tf_expr_group0)-1)/(length(tf_expr_group0))),4)
    sd_gene_expr_group0 = round(sd(gene_expr_group0)*sqrt((length(gene_expr_group0)-1)/(length(gene_expr_group0))),4)

    p0 = cor.test(tf_expr_group0 , gene_expr_group0 )
    r0 = signif(p0$estimate,4)
    p0 = signif(p0$p.value,3)
    
    # group 1
    mean_gene_expr_group1 = round(mean(gene_expr_group1),4)
    mean_tf_expr_group1 = round(mean(tf_expr_group1),4)
    sd_tf_expr_group1 = round(sd(tf_expr_group1)*sqrt((length(tf_expr_group1)-1)/(length(tf_expr_group1))),4)
    sd_gene_expr_group1 = round(sd(gene_expr_group1)*sqrt((length(gene_expr_group1)-1)/(length(gene_expr_group1))),4)

    p1 = cor.test(tf_expr_group1 , gene_expr_group1 )
    r1 = signif(p1$estimate)
    p1 = signif(p1$p.value,3)
    
    # group 2
    mean_gene_expr_group2 = round(mean(gene_expr_group2),4)
    mean_tf_expr_group2 = round(mean(tf_expr_group2),4)
    sd_tf_expr_group2 = round(sd(tf_expr_group2)*sqrt((length(tf_expr_group2)-1)/(length(tf_expr_group2))),4)
    sd_gene_expr_group2 = round(sd(gene_expr_group2)*sqrt((length(gene_expr_group2)-1)/(length(gene_expr_group2))),4)

    p2 = cor.test(tf_expr_group2 , gene_expr_group2 )
    r2 = signif(p2$estimate)
    p2 = signif(p2$p.value,3)
    
    gene_anova = CalAnova(gene_expr_group0,gene_expr_group1,gene_expr_group2)
    tf_anova = CalAnova(tf_expr_group0,tf_expr_group1,tf_expr_group2)
    
    mean_tf_expr = paste0(mean_tf_expr_group0,'/',mean_tf_expr_group1,'/',mean_tf_expr_group2)
    mean_gene_expr = paste0(mean_gene_expr_group0,'/',mean_gene_expr_group1,'/',mean_gene_expr_group2)
    
    std_gene_expr = paste0(sd_gene_expr_group0,'/',sd_gene_expr_group1,'/',sd_gene_expr_group2 )
    std_tf_expr = paste0(sd_tf_expr_group0,'/',sd_tf_expr_group1,'/',sd_tf_expr_group2 )
    
    counts = paste0(length(group_0),'/',length(group_1),'/',length(group_2))
    pccs = paste0(p0,'/',p1,'/',p2)
    cors  = paste0(r0,'/',r1,'/',r2)
    
    result = rbind(result,c(snp_gene[i,1],as.character(snp_gene[i,2]),snp_gene[i,3],mean_gene_expr,std_gene_expr,gene_anova,mean_tf_expr, std_tf_expr,tf_anova,counts,cors,pccs))
  }
  result = as.data.frame(result)
  colnames(result) =c('enhancer','snpid','tfgene','mean_gene_expr','std_gene_expr','gene_anova','mean_tf_expr','std_tf_expr','tf_anova','Counts','PCCS-cor','PCCS-pvalue')
  result
}

CalAnova = function( expr_0, expr_1, expr_2) {
  
  diff = 0
  exprsnp = NULL
  expr_factor = NULL
  if(length(expr_0)>0) {
    exprsnp =append(exprsnp,expr_0)
    expr_factor = append(expr_factor,c(rep("AA",length(expr_0))))
    diff = diff+1
  }
  
  if(length(expr_1)>0) {
    exprsnp =append(exprsnp,expr_1)
    expr_factor = append(expr_factor,c(rep("AG",length(expr_1))))
    diff = diff+1
  }
  
  if(length(expr_2)>0) {
    exprsnp =append(exprsnp,expr_2)
    expr_factor = append(expr_factor,c(rep("GG",length(expr_2))))
    diff = diff+1
  }
  
  
  exprsnp<-as.numeric(exprsnp)
  exprsnp_data.0<-data.frame(exprsnp,SNP=factor(expr_factor))
  if(diff>1) {
    p_anova.1 = aov(exprsnp~SNP,data=exprsnp_data.0)
    p_anova.1 = signif(summary.aov(p_anova.1)[[1]]$'Pr(>F)'[[1]],4)
  } else {
    p_anova.1 = NULL
  } 
  
  p_anova.1
}

#################################################################  chow test function
library(methods)

pan_chowtest = function(infile,gfile,cfile,outfile,covfile = '', pcut = 0.05) {
  
  logger = function(..., level = 'INFO') {
    cat(paste0(level, ': ', paste(...), '\n'), file = stderr())
  }
  
  # allow ifelse to return NULL
  ifelse = function(condition, true, false) {
    if(condition) return(true)
    return(false)
  }
  
  read.table.inopts = function(infile, inopts, dup = NULL, try = FALSE) {
    inopts.default = function(key, default) list.get(inopts, key, default, check.names = TRUE)
    optrnames = inopts.default('rnames', TRUE)
    #optrnames = ifelse('rnames' %in% opts, ifelse(inopts$rnames, 1, NULL), 1)
    params = list(
      infile,
      #header      = ifelse('cnames' %in% opts, inopts$cnames, T),
      header      = inopts.default('cnames', TRUE),
      row.names   = ifelse(!is.null(dup), NULL, ifelse(optrnames, 1, NULL)),
      sep         = inopts.default('delimit', "\t"),
      check.names = F,
      quote       = inopts.default('quote', ""),
      skip        = inopts.default('skip', 0)
    )
    if (!try) {
      d = do.call(read.table, params)
    } else {
      d = tryCatch({do.call(read.table, params)}, error = function(e) {
        logger('Error encountered while read file:', infile, level = 'WARNING')
        logger(e, level = 'WARNING')
        return (NULL)
      })
    }
    if (is.null(dup) || !optrnames || is.null(d)) {
      #return (d)
    } else {
      rnames = as.character(as.vector(d[,1,drop = T]))
      if (dup == 'drop') {
        rindex = !duplicated(rnames)
        d = d[rindex, -1, drop = FALSE]
        rownames(d) = rnames[rindex]
      } else if (dup == 'mean') {
        d = aggregate(d[,-1,drop=FALSE], by = list(rnames), FUN = mean)
        rownames(d) = as.character(d[,1,drop=TRUE])
        d = d[,-1,drop=FALSE]
      } else  { # keep
        rownames(d) = make.unique(rnames)
        d = d[,-1,drop = FALSE]
      }
    }
    d
  }
  
  # format data.frame to output
  pretty.numbers = function(df, formats) {
    # remember set stringsAsFactors as FALSE for the dataframe!!
    if (nrow(df) == 0) {
      return(df)
    }
    allCols      = colnames(df)
    formatedCols = c()
    for (fcols in names(formats)) {
      if (fcols == '.') { # must be last element of formats
        cols = which(!allCols %in% formatedCols)
      } else {
        cols = unlist(strsplit(fcols, '..', fixed = T))
        formatedCols = c(formatedCols, cols)
      }
      cols = intersect(cols, allCols)
      df[, cols] = sprintf(formats[[fcols]], as.numeric(unlist(df[, cols])))
    }
    df
  }
  
  # format data.frame to output
  pretty.numbers2 = function(df, ...) {
    formats = list(...)
    options(stringsAsFactors = FALSE)
    df = as.data.frame(df)
    if (nrow(df) == 0)
      return(df)
    
    allCols      = colnames(df)
    if (is.null(allCols))
      allCols = 1:ncol(df)
    formatedCols = c()
    for (fcols in names(formats)) {
      if (fcols == '.') { # must be last element of formats
        if (length(formatedCols) == 0) {
          cols = allCols
        } else {
          cols = which(!allCols %in% formatedCols)
        }
      } else {
        cols = unlist(strsplit(fcols, '..', fixed = T))
        formatedCols = c(formatedCols, cols)
      }
      cols = intersect(cols, allCols)
      df[, cols] = sprintf(formats[[fcols]], as.numeric(unlist(df[, cols])))
    }
    df
  }
  
  is.installed = function(pkg) {
    is.element(pkg, installed.packages()[,1])
  }
  
  bQuote = function(s) {
    if (startsWith(s, '`') && endsWith(s, '`')) {
      return (s)
    } else {
      paste0('`', s, '`')
    }
  }
  
  is.true = function(x, collapse = 'all') {
    if (is.null(x)) return (FALSE)
    if (length(x) == 0) return (FALSE)
    if (length(x) == 1) {
      if (is.na(x)) return (FALSE)
      if (is.list(x)) return (TRUE)
      if (is.character(x)) return (nchar(x) > 0)
      tryCatch({
        x = as.logical(x)
      }, error = function(e){
        x <<- TRUE
      })
      if (is.na(x)) return (TRUE)
      return (x)
    } else if (collapse == 'all') {
      for (i in x) {
        if (!is.true(i, 'any')) return (FALSE)
      }
      return (TRUE)
    } else {
      for (i in x) {
        if (is.true(i, 'any')) return (TRUE)
      }
      return (FALSE)
    }
  }
  
  is.false = function(x, collapse = 'all') {
    !is.true(x, ifelse(collapse == 'all', 'any', 'all'))
  }
  
  list.get = function(l, key, default = NULL, check.names = FALSE) {
    # get the value of a key in list with default
    # @params:
    #	`l`: The list
    #	`key`: The key
    #	`default`: The default value. Default: `NULL`
    #	`check.names`: Check whetheer the name exists, even with value `NULL`. Default: `FALSE`
    #		- `list.get(list(a = NULL), 'a', default = 1, check.names = TRUE) == NULL`
    #		- `list.get(list(a = NULL), 'a', default = 1, check.names = FALSE) == 1`
    if (!check.names) {
      ifelse(is.null(l[[key]]), default, l[[key]])
    } else {
      ns = names(l)
      if (key %in% ns)
        return (l[[key]])
      return (default)
    }
  }
  options(stringsAsFactors = F)
  outdir   = ""
  dofdr    = TRUE
  plotchow = FALSE
  devpars  = ''#Box(res = 300, width = 2000, height = 2000)
  ggs      = ''#box()
  inopts   =  ''#box(cnames = True, rnames = True)
  # covfile  = ''
  
  if (plotchow) {
    {{rimport}}('plot.r')
  }
  
  if (dofdr == T) dofdr = 'BH'
  
  chow.test = function(formula, group, data, covdata = NULL, ...) {
    fmvars  = all.vars(as.formula(formula))
    vars    = colnames(data)
    if (length(fmvars) == 2 && fmvars[2] == '.') {
      fmvars = c(fmvars[1], vars[vars!=fmvars[1] & vars!=group])
    }
    formula = sprintf(
      '%s ~ %s',
      bQuote(fmvars[1]),
      paste(sapply(fmvars[2:length(fmvars)], bQuote), collapse = '+')
    )
    covs = NULL
    if (is.null(covdata)) {
      pooledfm = as.formula(formula)
    } else {
      covdata = covdata[rownames(data),,drop = FALSE]
      covs    = colnames(covdata)
      data    = cbind(data, covdata)
      rm(covdata)
      pooledfm = as.formula(paste(formula, '+', paste(sapply(covs, bQuote), collapse = '+')))
      fmvars   = c(fmvars, covs)
    }
    
    if (sum(complete.cases(data[,fmvars])) < 2) {
      pooled_lm = NULL
    } else {
      pooled_lm = lm(pooledfm, data = data, ...)
    }
    #coeff     = as.list(pooled_lm$coefficients)
    groups    = levels(as.factor(data[, group]))
    group_lms = sapply(groups, function(g) {
      if (is.null(covdata)) {
        subfm  = as.formula(formula)
      } else {
        subfm = as.formula(paste(
          formula, '+', 
          #paste(sapply(covs, function(x) paste0('offset(', coeff[[x]], '*', bQuote(x), ')')), collapse = '+')
          paste(sapply(covs, bQuote), collapse = '+')
        ))
      }
      sublmdata = data[data[,group] == g, , drop = FALSE]
      if (sum(complete.cases(sublmdata[,fmvars])) < 2) {
        NULL
      } else {
        list(lm(subfm, data = sublmdata, ...))
      }
    })
    
    pooled.ssr = ifelse(is.null(pooled_lm), NA, sum(pooled_lm$residuals ^ 2))
    subssr     = ifelse(is.false(group_lms, 'any'), NA, sum(sapply(group_lms, function(x) sum(x$residuals ^ 2))))
    ngroups    = length(groups)
    K          = length(fmvars) + length(covs)
    J          = (ngroups - 1) * K
    DF         = nrow(data) - ngroups * K
    FS         = (pooled.ssr - subssr) * DF / subssr / J
    list(
      pooled.lm  = pooled_lm,
      group.lms  = group_lms,
      Fstat      = FS,
      group      = group,
      pooled.ssr = pooled.ssr,
      group.ssr  = subssr,
      Pval       = pf(FS, J, DF, lower.tail = FALSE)
    )
  }
  
  plot.chow = function(chow, plotfile, ggs, devpars) {
    cols     = all.vars(chow$pooled.lm$terms)[1:2]
    plotdata = do.call(rbind, lapply(names(chow$group.lms), function(m) data.frame(chow$group.lms[[m]]$model[, cols, drop = FALSE], group = m)))
    colnames(plotdata)[3] = chow$group
    if (!is.null(ggs$scale_color_discrete)) {
      ggs$scale_color_discrete$name = ifelse(
        is.function(ggs$scale_color_discrete$name), 
        ggs$scale_color_discrete$name(chow$group),
        chow$group
      )
      ggs$scale_color_discrete$labels = sapply(names(chow$group.lms), function(x) {
        coeff = as.list(chow$group.lms[[x]]$coefficients)
        bquote(.(x): beta[.(cols[2])]==.(round(coeff[[cols[2]]], 3)) ~ "," ~ epsilon == .(round(coeff[['(Intercept)']], 3)))
      })
    }
    
    if (!is.null(ggs$scale_shape_discrete)) {
      ggs$scale_shape_discrete$name = ifelse(
        is.function(ggs$scale_shape_discrete$name), 
        ggs$scale_shape_discrete$name(chow$group),
        chow$group
      )
      ggs$scale_shape_discrete$labels = sapply(names(chow$group.lms), function(x) {
        coeff = as.list(chow$group.lms[[x]]$coefficients)
        bquote(.(x): beta[.(cols[2])]==.(round(coeff[[cols[2]]], 3)) ~ "," ~ epsilon == .(round(coeff[['(Intercept)']], 3)))
      })
    }
    
    if (is.null(ggs$scale_color_discrete) && is.null(ggs$scale_shape_discrete)) {
      ggs$scale_color_discrete = list(
        name = chow$group,
        labels = sapply(names(chow$group.lms), function(x) {
          coeff = as.list(chow$group.lms[[x]]$coefficients)
          bquote(.(x): beta[.(cols[2])]==.(round(coeff[[cols[2]]], 3)) ~ "," ~ epsilon == .(round(coeff[['(Intercept)']], 3)))
        })
      )
    }
    
    if (is.null(ggs$scale_color_discrete)) ggs$scale_color_discrete = ggs$scale_shape_discrete
    if (is.null(ggs$scale_shape_discrete)) ggs$scale_shape_discrete = ggs$scale_color_discrete
    
    plot.points(
      plotdata, 
      plotfile, 
      x = 2, y = 1, 
      params = list(aes_string(color = chow$group, shape = chow$group)), 
      ggs = c(ggs, list(
        geom_smooth = list(aes_string(color = chow$group), method = "lm", se = FALSE)
      ))
    )
  }
  
  formatlm = function(m) {
    if (class(m) == 'lm') {
      coeff = as.list(m$coefficients)
      vars = all.vars(m$terms)
      terms = unlist(sapply(c(vars[2:length(vars)], '(Intercept)', 'N'), function(x) {
        ce = list.get(coeff, x, list.get(coeff, bQuote(x)))
        if (x == 'N') {
          paste0('N=', nrow(m$model))
        } else if (is.null(ce)) {
          NULL
        } else {
          l = ifelse(x == '(Intercept)', '_', x)
          paste0(l, '=', round(ce, 3))
        }
      }))
      paste(terms[!is.null(terms)], collapse = ', ')
    } else {
      paste(sapply(names(m), function(x) {
        paste0(x, ': ', formatlm(m[[x]]))
      }), collapse = ' // ')
    }
  }
  
  results = data.frame(
    Case   = character(),
    Pooled = character(),
    Groups = character(),
    SSR    = double(),
    SumSSR = double(),
    Fstat  = double(),
    Pval   = double()
  )
  
  indata = read.table.inopts(infile, inopts, try = TRUE)
  if (is.null(indata)) {
    write.table(results, outfile, col.names = T, row.names = F, sep = "\t", quote = F)
    quit(save = "no")
  }
  
  #     X1  X2  X3  X4 ... Y
  # G1  1   2   1   4  ... 9
  # G2  2   3   1   1  ... 3
  # ... ...
  # Gm  3   9   1   7  ... 8
  #K      = ncol(indata)
  covdata = NULL
  covs    = NULL
  if (covfile != "") {
    covdata = read.table(covfile, header = T, row.names = 1, check.names = F)
    #indata  = cbind(covdata[rownames(indata),,drop = F], indata)
    covs = colnames(covs)
  }
  gdata  = read.table.inopts(gfile, list(cnames = TRUE, rnames = TRUE))
  # 	Case1	Case2
  # G1	Group1	Group1
  # G2	Group1	NA
  # G3	Group2	Group1
  # ... ...
  # Gm	Group2	Group2
  cases  = colnames(gdata)
  fmulas = data.frame(x = rep(paste(bQuote(colnames(indata)[ncol(indata)]), '~ .'), length(cases)))
  rownames(fmulas) = cases
  if (!is.null(cfile) && cfile != "") {
    fmulas = read.table.inopts(cfile, list(cnames = FALSE, rnames = TRUE))
    cases  = rownames(fmulas)
  }
  
  for (case in cases) {
    logger('Handling case: ', case, '...')
    fmula  = fmulas[case,,drop = TRUE]
    groups = gdata[!is.na(gdata[,case]),case,drop = FALSE]
    data   = cbind(indata[rownames(groups),, drop = FALSE], group = groups)
    colnames(data)[ncol(data)] = case
    ct = chow.test(fmula, case, data, covdata = covdata)
    if (dofdr == FALSE && (is.na(ct$Pval) || ct$Pval >= pcut)) {
      next
    }
    results = rbind(results, list(
      Case   = case,
      Pooled = formatlm(ct$pooled.lm),
      Groups = formatlm(ct$group.lms),
      SSR    = ct$pooled.ssr,
      SumSSR = ct$group.ssr,
      Fstat  = ct$Fstat,
      Pval   = ct$Pval
    ))
    # doplot
    if (plotchow && ct$Pval < pcut) {
      plot.chow(ct, file.path(outdir, paste0(case, '.png')), ggs, devpars)
    }
  }
  
  if (dofdr != F) {
    results = cbind(results, Qval = p.adjust(results$Pval, method = dofdr))
  } 
  write.table(pretty.numbers(results, list(
    SSR..SumSSR..Fstat = '%.3f',
    Pval..Qval = '%.3E'
  )), outfile, col.names = T, row.names = F, sep = "\t", quote = F)
}

### main ###
if(DEBUG){
  tcga_path = '/research/labs/pharmacology/junwenwang/data/yanxi/tcga_test/'
  fgeno = "/research/labs/pharmacology/junwenwang/data/yanxi/tcga_test/snp_mtx.txt"
  fexpr = "/research/labs/pharmacology/junwenwang/data/yanxi/tcga_test/expr_mtx.txt"
  fgene_loc = "/research/labs/pharmacology/junwenwang/data/yanxi/reference/TCGA_gene_loc_strand_hg38.txt"
  ftf_info = "/research/labs/pharmacology/junwenwang/data/genecard/all/info/TP73.txt"
  fcov = "/research/labs/pharmacology/junwenwang/data/yanxi/tcga_test/cov.txt"
  genename = "TP73"
  
  ceQTL(genename, fgeno, fexpr, fgene_loc, ftf_info, fcov)
} else{
  args <- commandArgs(trailingOnly = TRUE)
  #tcga_path = '/research/labs/pharmacology/junwenwang/data/tcga/GBM_genecard/'
  tcga_path = as.character(args[1])
  #setwd(tcga_path)
  snp = read.csv(paste0(tcga_path,as.character(args[2])),header = T,sep="\t")   # race
  expr = read.csv(paste0(tcga_path,as.character(args[3])),header = T,sep="\t")
  genotype = read.table(paste0(tcga_path,as.character(args[4])),header = T,sep = '\t')
  gene_loc = read.table(paste0(tcga_path,as.character(args[5])),header = T,sep='\t')
  snp_loc = read.table(paste0(tcga_path,as.character(args[6])),header = T,sep='\t')

  #gene_list = read.table('/research/labs/pharmacology/junwenwang/data/genecard/gpcr/info/XCR1.txt',header = T,sep = '\t')
  #gene_list = read.table('./info/ETV1.txt',header = T,sep = '\t')
  source_path =  as.character(args[7])
  gene_list = read.table(paste0(source_path,as.character(args[8])),header = T,sep = '\t')#  MSTN_gene.txt \ EMP2.txt\THY1.txt\TF.txt
  gene_list = gene_list[,c(1,5,4,7)]
  colnames(gene_list) = c("Enhancer_id","Loc","tss_dis","TFs" )
  gene_list$TFs = as.character(gene_list$TFs)
}

