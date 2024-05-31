library(ggplot2)
library(tidyverse)
library(viridis)
library(ggpubr)
library(readxl)

tcga_cor_plot <- function(filename, sep, ref, alt, geneX, geneY, snp, code, out, model = "default"){
  data <- read.csv(filename, sep = sep)
  
  data.wild <- data[which(data$genotype == paste(ref, ref, sep = "")), ]
  data.het <- data[which(data$genotype == paste(ref, alt, sep = "")), ]
  data.homo <- data[which(data$genotype == paste(alt, alt, sep = "")), ]
  data.dom <- data[which(data$genotype != paste(ref, ref, sep = "")), ]
  data.rec <- data[which(data$genotype != paste(alt, alt, sep = "")), ]
  
  # linear regression plot with pearson and p-value
  if(out == ""){
    out = getwd()
  }
  if (!dir.exists(file.path(getwd(), out))){
    dir.create(file.path(getwd(), out))
  }
  if(model == "default"){
    pdf(file.path(out, paste(code, '-', geneX, '-', snp, '-', geneY, '_default.pdf', sep = "")))
    # jpeg(paste("correlation_plot_", geneX, geneY, ".jpg", sep = ""))
    
    # plt.all <- ggplot(data, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
    #   geom_point()+
    #   geom_smooth(method=lm) +
    #   stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson")
    
    plt.wild <- ggplot(data.wild, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
      ggtitle(paste(paste(ref, ref, sep = ""), "(n=", nrow(data.wild), ")", sep = "")) +
      geom_point()+
      geom_smooth(method=lm) +
      theme(plot.title = element_text(hjust = 0.5)) +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), fontface = "italic"), method = "pearson") +
      labs(x = paste0(geneX, "(RSEM)"), y = paste0(geneY, "(RSEM)"))
    
    plt.het <- ggplot(data.het, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
      ggtitle(paste(paste(ref, alt, sep = ""), "(n=", nrow(data.het), ")", sep = "")) +
      geom_point()+
      geom_smooth(method=lm) +
      theme(plot.title = element_text(hjust = 0.5)) +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), fontface = "italic"), method = "pearson") +
      labs(x = paste0(geneX, "(RSEM)"), y = paste0(geneY, "(RSEM)"))
    
    plt.homo <- ggplot(data.homo, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
      ggtitle(paste(paste(alt, alt, sep = ""), "(n=", nrow(data.homo), ")", sep = "")) +
      geom_point()+
      geom_smooth(method=lm) +
      theme(plot.title = element_text(hjust = 0.5)) +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), fontface = "italic"), method = "pearson") +
      labs(x = paste0(geneX, "(RSEM)"), y = paste0(geneY, "(RSEM)"))
    
    plt <- ggarrange(plt.wild, plt.het, plt.homo,
                     # labels = c(paste(ref, ref, sep = ""), paste(ref, alt, sep = ""), paste(alt, alt, sep = "")),
                     # font.label = list(size = 10, face = "plain"),
                     ncol = 2, nrow = 2)
    print(plt)
    dev.off()
  } else if(model == "dominant"){
    pdf(paste(out, "/", "plots", "/", code, "-", geneX, "-", snp, "-", geneY, "_dominant.pdf", sep = ""))
    # jpeg(paste("correlation_plot_", geneX, geneY, ".jpg", sep = ""))
    
    # plt.all <- ggplot(data, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
    #   geom_point()+
    #   geom_smooth(method=lm) +
    #   stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson")
    
    plt.wild <- ggplot(data.wild, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
      ggtitle(paste(paste(ref, ref, sep = ""), "(n=", nrow(data.wild), ")", sep = "")) +
      geom_point()+
      geom_smooth(method=lm) +
      theme(plot.title = element_text(hjust = 0.5)) +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), fontface = "italic"), method = "pearson") +
      labs(x = paste0(geneX, "(RSEM)"), y = paste0(geneY, "(RSEM)"))
    
    plt.dom <- ggplot(data.dom, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
      ggtitle(paste(paste(ref, alt, "/", alt, alt, sep = ""), "(n=", nrow(data.dom), ")", sep = "")) +
      geom_point()+
      geom_smooth(method=lm) +
      theme(plot.title = element_text(hjust = 0.5)) +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), fontface = "italic"), method = "pearson") +
      labs(x = paste0(geneX, "(RSEM)"), y = paste0(geneY, "(RSEM)"))
    
    plt <- ggarrange(plt.wild, plt.dom,
                     # labels = c(paste(ref, ref, sep = ""), paste(ref, alt, sep = ""), paste(alt, alt, sep = "")),
                     # font.label = list(size = 10, face = "plain"),
                     ncol = 2, nrow = 2)
    print(plt)
    dev.off()
  }else if(model == "recessive"){
    pdf(paste(out, "/", "plots", "/", code, "-", geneX, "-", snp, "-", geneY, "_recessive.pdf", sep = ""))
    # jpeg(paste("correlation_plot_", geneX, geneY, ".jpg", sep = ""))
    
    # plt.all <- ggplot(data, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
    #   geom_point()+
    #   geom_smooth(method=lm) +
    #   stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson")
    
    plt.rec <- ggplot(data.rec, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
      ggtitle(paste(paste(ref, ref, "/", ref, alt, sep = ""), "(n=", nrow(data.rec), ")", sep = "")) +
      geom_point()+
      geom_smooth(method=lm) +
      theme(plot.title = element_text(hjust = 0.5)) +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), fontface = "italic"), method = "pearson") +
      labs(x = paste0(geneX, "(RSEM)"), y = paste0(geneY, "(RSEM)"))
    
    plt.homo <- ggplot(data.homo, aes(x=.data[[geneX]], y=.data[[geneY]])) + 
      ggtitle(paste(paste(alt, alt, sep = ""), "(n=", nrow(data.homo), ")", sep = "")) +
      geom_point()+
      geom_smooth(method=lm) +
      theme(plot.title = element_text(hjust = 0.5)) +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), fontface = "italic"), method = "pearson") +
      labs(x = paste0(geneX, "(RSEM)"), y = paste0(geneY, "(RSEM)"))
    
    plt <- ggarrange(plt.rec, plt.homo,
                     # labels = c(paste(ref, ref, sep = ""), paste(ref, alt, sep = ""), paste(alt, alt, sep = "")),
                     # font.label = list(size = 10, face = "plain"),
                     ncol = 2, nrow = 2)
    print(plt)
    dev.off()
  } else{
    print("model error!")
  }

  return(plt)
}

# # box plot
# pdf("box_plt_VDR.pdf")
# # jpeg("box_plt_VDR.jpg")
# box_plt_VDR <- ggboxplot(data, x='genotype', y='VDR', fill='genotype') + 
#   # geom_boxplot() + 
#   geom_jitter(color="black", size=0.4, alpha=0.9) +
#   stat_compare_means(method = "anova")
# box_plt_VDR
# dev.off()

### main ### 
# filename = "/Users/ychen745/Desktop/plt_test/snp_list.txt"
args = commandArgs(trailingOnly = TRUE)
filename = as.character(args[1])
input_folder = as.character(args[2])
out = as.character(args[3])
model = as.character(args[4])
if(substr(input_folder, nchar(input_folder), nchar(input_folder)) == "/"){
  input_folder = substr(input_folder, 1, nchar(input_folder)-1)
}
print(out)
if(substr(out, nchar(out), nchar(out)) == "/"){
  out = substr(out, 1, nchar(out)-1)
}
dir.create(file.path(out, "plots"), showWarnings = FALSE)
path_list = unlist(strsplit(input_folder, "/"))
folder = paste(path_list[1:length(path_list)-1], collapse = "/")
info_mtx = read.csv(filename, sep = "\t", colClasses = replicate(6, "character"))
delimiter = "\t"
models = c("default", "dominant", "recessive")
for(i in 1:nrow(info_mtx)){
  for(model in models){
    rsid = info_mtx$rsid[i]
    enhancer = info_mtx$enhancer[i]
    target = info_mtx$target[i]
    ref = info_mtx$ref[i]
    alt = info_mtx$alt[i]
    code = info_mtx$code[i]
    fname = paste(input_folder, "/", paste(c(code, target, rsid, enhancer), collapse = "-"), ".txt", sep = "")
    print(fname)
    # print(fname)
    tcga_cor_plot(fname, delimiter, ref, alt, enhancer, target, rsid, code, out, model)
  }
}
############

### debugging ###
# filename = "plt_20220308/01-SOAT1-chr1_179295262_rs10913708_A_G-ABL2.txt"
# ref = "G"
# alt = "A"
# enhancer = "ABL2"
# target = "SOAT1"
# rsid = "chr1_179295262_rs10913708_A_G"
# code = "01"
# out = "plt_result_20220308"
# 
# tcga_cor_plot(filename, "\t", ref, alt, enhancer, target, rsid, code, out, "recessive")
### debugging ###
# pdf("test.pdf")
# print(plt)
# dev.off()
