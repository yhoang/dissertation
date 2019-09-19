#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
rm(list = ls())
source("R/spitzer.r")


alphalist <- 0.9
checkCALC = "range.quad"
tissue ="BL"
iterations = 500
quad.cap = 5000
bin.min.count = 400
top = 10

folder = sprintf("table/",tissue,checkCALC,alphalist)
file.name = paste0(sprintf("%s_%s_signif_quadrants%s_it%s_table.csv",tissue,checkCALC,quad.cap,iterations))

tissue = c("BL")

cofactor = 0.1
plot.folder = folder

marker.comb = read.table(paste0(folder,file.name),sep="\t",header=T)[1:top,2]
#marker.comb = read.table(paste0(folder,file.name),sep=",",header=T)[1:10,2]


for ( type in tissue) {
  
  print(type)
  if (type == "BL") {
    tissue1="blood"
    file.name1 = sprintf("Rdata/d3_%s_triploT_quadrants_%s_cof%s_minbin%s_keep_bio_cut0_trimmed.rds",
                         tissue1,checkCALC,cofactor,bin.min.count)
    df1 = readRDS(file = file.name1)
    metadata1 = read.csv(file=sprintf("meta/Spitzer_%s.csv",tissue1,tissue1))
    df1.condition = metadata1$condition2[which(as.character(metadata1$sample_id) %in% rownames(df1) )]
    df1= cbind(df1.condition,df1)
    names(df1)[1] = "condition"
    dim(df1)
    
    
    ### remove untr_d3_Bl3
    df1 = df1[-13,]
    
    df.total = df1
    df.total = as.matrix(df.total)
    ## convert NaN/+-Inf to -1
    df.total[is.nan(df.total) | is.infinite(df.total)] <- -1
    
  } else {
    df.total = readRDS(sprintf("/scratch/drfz/Spitzer2017/Rdata_cof/d3_%s_triploT_quadrants_%s_cof%s_bio_cut0_trimmed_NaN-1_-Q1.rds",type,checkCALC,cofactor))
    #df.total = readRDS(sprintf("/scratch/drfz/Spitzer2017/Rdata_cof/d3_%s_triploT_quadrants_%s_cof%s_bio_cut0_NaN0.rds",type,checkCALC,cofactor))
    dim(df.total)
  }
  ## remove columns with "Q1" (10962 variables)
  if (any(grepl("Q1",colnames(df.total)))) df.total = df.total[,-grep("Q1",colnames(df.total))]

  print(dim(df.total))
  #[1]    13 2523
  
  sample.size = ncol(df.total)-1
  
  
  plot.file = sprintf("%s/%s_%s_boxplots_top%s_%s.pdf",plot.folder,type,checkCALC,top,sample.size)
  
  df.total = as.data.frame(df.total)
  quad.idx = which(colnames(df.total) %in% marker.comb)
  df.predict = df.total[order(rownames(df.total)),c(1,quad.idx)]
  
  ### manipulate dataframe
  dat.m = melt(df.predict, id.vars="condition")
  
  ### save boxplot as pdf
  plot.pred = ggboxplot(dat.m, x="condition",y="value",ylim = c(0,6),
                        color = "condition", palette=c("#b2182b","#006837"), add = "jitter",
                        facet.by = "variable",ncol = 4,
                        xlab="",width=0.3
                        ) +
    theme(text = element_text(size=12),
          axis.title.x.top = element_text(colour="grey20",
                                          size=10,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.x.bottom = element_text(colour="black",size=12),
          plot.margin = unit(c(0.5,1,0.3,1), "cm")) +
    stat_compare_means(method="t.test",label = "p.signif", label.x = 1.5, label.y = 0)
  ggsave(plot.file,plot.pred,width=15,height=25)
}
printf("Written in %s",plot.file)
