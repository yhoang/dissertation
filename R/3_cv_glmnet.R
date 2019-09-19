#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
rm(list = ls())
source("R/spitzer.r")


set.cluster = 1
alphalist <- 0.9
checkCALC = "range.quad"
tissue ="BL"
iterations = 500
quad.cap = 5000
cofactor = 0.1
bin.min.count = 400
folder = sprintf("table/",tissue,checkCALC,alphalist)

### helpful functions
printf <- function(...) invisible(print(sprintf(...)))
is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}

tissue1="blood"
df.total = readRDS(file = sprintf("Rdata/d3_%s_triploT_quadrants_%s_cof%s_minbin%s_keep_bio_cut0_trimmed.rds",
                                  tissue1,checkCALC,cofactor,bin.min.count))

metadata1 = read.csv(file=sprintf("meta/Spitzer_blood.csv",tissue1,tissue1))
df1.condition = metadata1$condition2[which(as.character(metadata1$sample_id) %in% rownames(df.total) )]
df.total= cbind(df1.condition,df.total)
names(df.total)[1] = "condition"

### remove untr_d3_Bl3
df.total=df.total[-13,]
dim(df.total)
#[1]   13 2523

df.total = as.matrix(df.total)
### convert NaN/+-Inf to -1
df.total[is.nan(df.total) | is.infinite(df.total)] <- -1
### convert NAs to sample group mean
for ( i in 2:ncol(df.total)) {
  NA.idx = which(is.na(df.total[,i]))
  
  for (j in NA.idx) {
    tmp = df.total[which(df.total[j,1]==df.total[,1]),i]
    tmp = tmp[-which(is.na(tmp))]
    tmp = round(sample(seq(mean(tmp)-sd(tmp),mean(tmp)+sd(tmp),by=0.01),1),2)
    
    df.total[j,i] = tmp
  }
}

sample.size = ncol(df.total)-1

### column index of predicted variable in dataset
typeColNum = 1
it_rmse = 0
it_run = 0
conf_interval=0.95

acc_total = vector()
padj_total = 0

seed.vec = sample(10^5)

#cl <- makeCluster(set.cluster, manual=TRUE, outfile="")
#cl <- makeCluster(set.cluster)
# cl<-parallel::makeCluster(set.cluster, verbose = TRUE)
# registerDoParallel(cl)

timeSTART = Sys.time()
ptm <- proc.time()
printf("##### Start at %s.",timeSTART)
printf("run %s::%s::samplesize=%s",tissue,checkCALC,sample.size)
it.total = it.run = 0
while (it.run < iterations) {
  
  it.total = it.total + 1  
  
  ### split into training and test sets
  set.seed(seed.vec[it.total]); train.idx = sample(nrow(df.total),ceiling(nrow(df.total)*0.8))
  trainset = df.total[train.idx,]
  set.seed(seed.vec[it.total]); set.foldid = sample(rep(seq((1/3)*nrow(trainset)),length=nrow(trainset)))
  testset = df.total[-train.idx,]
  
  if (
    all(trainset[which(set.foldid==1),typeColNum]==0) |
    all(trainset[which(set.foldid==2),typeColNum]==0) |
    all(trainset[which(set.foldid==3),typeColNum]==0) |
    all(trainset[which(set.foldid==1),typeColNum]==1) |
    all(trainset[which(set.foldid==2),typeColNum]==1) |
    all(trainset[which(set.foldid==3),typeColNum]==1) |
    all(testset[,typeColNum]==0) |
    all(testset[,typeColNum]==1) 
  ) next;
  
  it.run = it.run + 1
  
  options(na.action = "na.pass")
  ### no need to set nfolds if foldid is provided since observations are already distributed in folds
  cv.out <- cv.glmnet(x = as.matrix(trainset[,-typeColNum]), y = trainset[,typeColNum],
              alpha=alphalist, family="binomial",
              lambda.min.ratio=.0005,
              type.measure="deviance",
              foldid = set.foldid
              # ,parallel = TRUE
  )
  
  # plot result
  # par(oma=c(1,1,1,1),mar=c(1,1,1,1))
  # plot(cv.out)

  # min value of lambda
  lambda_min <- cv.out$lambda.min
  # best value of lambda
  lambda_1se <- cv.out$lambda.1se
  
  ### prediction
  lasso_class <- predict(cv.out,
                         newx = as.matrix(testset[,-typeColNum]),
                         s=lambda_min,
                         na.action = na.pass,
                         type="class"
                 )
  lasso_prob <- predict(cv.out,newx = as.matrix(testset[,-typeColNum]),s=lambda_1se,type="response")
  lasso_train <- predict(cv.out,newx = as.matrix(trainset[,-typeColNum]),s=lambda_1se,type="response")
  
  ### accuracy
  acc = mean(lasso_class==testset[,typeColNum])
  acc_total = c(acc_total,acc)
  
  ### root mean squared error RMSE
  # in sample RMSE
  rmse_train = sqrt(sum(lasso_train - mean(trainset[,typeColNum]))^2/nrow(trainset))
  # out of sample RMSE
  rmse_test = sqrt(sum(lasso_prob - mean(testset[,typeColNum]))^2/nrow(testset))
  # diff
  rmse_diff = rmse_test-rmse_train
  
  
  ### regression coefficients
  coeff = coef(cv.out,s=lambda_1se)

  ### get coefficients who has impact
  coeff.idx = which(abs(coeff)>0)
  
  coeff.pred = coeff[coeff.idx,]
  coeff.pred = coeff.pred[order(abs(coeff.pred),decreasing=T)]
  
    if ((rmse_diff)<0.05 & length(coeff.pred)>1) {
      it_rmse = it_rmse+1
      if (length(coeff.pred) > quad.cap) coeff.pred = coeff.pred[1:(quad.cap+1)]
      
      ### look at differentiating variables (here = quadrants)
      intercept.idx = grep("Intercept",names(coeff.pred))
      vars.pred = names(coeff.pred)[-intercept.idx]
      vars.pred = gsub("`","",vars.pred)
      
      df.total = as.data.frame(df.total)
      df.predict = df.total[order(rownames(df.total)),c(1,which(names(df.total)%in%vars.pred))]
      
      ### manipulate dataframe
      dat.m = melt(df.predict, id.vars="condition")
      
      ### unpaired two samples t-test (independent samples)
      ### calculate t.test for each variable
      comp = compare_means(value ~ condition, data = dat.m, group.by = "variable",  
                           method="t.test", p.adjust.method = "BH"
                            )
      
      ### get idx where adjusted p-value is <0.05
      comp.adj.idx = which(comp$p.adj<0.05)
      padj_total = padj_total+length(comp.adj.idx)
      
      ### get quadrants and coeffs with p.adj<0.05
      coeff.pred.adj = coeff.pred[which(names(coeff.pred) %in% comp$variable[comp.adj.idx] )]
      
      ### write info to signif_quadrants.txt
      # create file and header
      write(paste("tissue","calculation_meth","sample_size","prediction_RMSE","seed","quadrants_padj",sep="\t"),
       file=sprintf("%ssignif_quadrants.txt",folder),append=TRUE)
      write(paste(tissue,checkCALC,"quads",round(rmse_test,4),seed.vec[it.total],
                  paste(names(coeff.pred.adj),collapse="\t"),
                  sep="\t"),file=sprintf("%ssignif_quadrants%s_it%s.txt",folder,quad.cap,iterations),append=TRUE)
      write(paste(tissue,checkCALC,"coeff",round(rmse_test,4),seed.vec[it.total],
                  paste(coeff.pred.adj,collapse="\t"),
                  sep="\t"),file=sprintf("%ssignif_quadrants%s_it%s.txt",folder,quad.cap,iterations),append=TRUE)
      
      
      ### get protein names    
      protein.names = readRDS("/scratch/drfz/Spitzer2017/protein_names")
      
      ### count protein in protein combination of quadrants and calculation method
      var.list = sapply(protein.names, grep, comp$variable[comp.adj.idx])#,value=T)
      #var.list[which( lengths(var.list) != 0)]
      #lengths(var.list[which( lengths(var.list) != 0)])
      #unlist(sapply(protein.names, grepl, comp$variable[comp.adj.idx]))
      var.count = lengths(var.list)
      
      ### write only counts of proteins from significant quadrants
      # create file and header
      write(paste("tissue","calculation_meth","sample_size","prediction_RMSE","seed",
                 paste(protein.names,collapse="\t"),sep="\t"),
                 file=sprintf("%svars_of_signif_quadrants.txt",folder),append=TRUE)
      write(paste(tissue,checkCALC,sample.size,round(rmse_test,2),seed.vec[it.total],paste(var.count,collapse="\t"),sep="\t"),
            file=sprintf("%svars_of_signif_quadrants%s_it%s.txt",folder,quad.cap,iterations),append=TRUE)
    } else {
      write(paste(tissue,checkCALC,sample.size,round(rmse_test,2),seed.vec[it.total],sep="\t"),
            file=sprintf("%ssignif_quadrants%s_it%s.txt",folder,quad.cap,iterations),append=TRUE)
      write(paste(tissue,checkCALC,sample.size,round(rmse_test,2),seed.vec[it.total],sep="\t"),
            file=sprintf("%svars_of_signif_quadrants%s_it%s.txt",folder,quad.cap,iterations),append=TRUE)
    }
  
  
  
  it_run = it_run + 1
}

printf("##### Ran %scv_glmnet.r with %s iterations and cap(quad)=%s. End at %s.",folder,it.total,quad.cap,Sys.time())
print(difftime(Sys.time(),timeSTART))

printf("Ready for glmnet alpha%s::%s::%s::%s.",alphalist,tissue,checkCALC,it.total)

printf("acc_total=%s, padj_total=%s",table(acc_total),padj_total)
printf("it_run=%s, it_rmse<0.5=%s",it_run,it_rmse)

print(proc.time() - ptm)