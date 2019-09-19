#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
rm(list = ls())

############################## Libraries
### for data access in PRI-base
library(RSQLite)
### for mds plotting
library(limma)
### for elastic-net regularization and cross-validation
library(glmnet)
### for data manipulation
library(reshape2)
library(dplyr)
### for parallel use (optional)
library(foreach)
library(doParallel)

### for plotting
library(ggplot2)
# https://www.r-bloggers.com/add-p-values-and-significance-levels-to-ggplots/
# devtools::install_github("kassambara/ggpubr")
library(ggpubr)
##############################


### create new environment
#rm(list = ls())
fcs = new.env()
fcs$parent.env=ls()

### Helpful functions
printf <- function(...) invisible(print(sprintf(...)))
is.even <- function(x) x %% 2 == 0
is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}

############################## Data import
###### connects to database fname
# @fname  file name of data base
fcs$connectDb<- function(fname) {
  this=fcs
  this$conn <- dbConnect(SQLite(), dbname = fname)
  print(paste("Database opened:",fname))
}

###### disconnect to database this$conn
fcs$disconnectDb <- function (dname) {
  this=fcs
  dbDisconnect(this$conn)
  print(paste("Database closed:",file.path(dname)))
}

##### Get meta table from database
##### @param
# table         project name
fcs$getDFtable <- function (table) {
  this=fcs
  table.df=dbGetQuery(this$conn,paste("SELECT * FROM ", table))
  
  return(table.df)
}

##### Get marker names for selected file
##### @param
# staintable    stain vector (optional)
# index         file index (optional)   
fcs$getVariables <- function (
  staintable=NA,
  table=NA,
  index=NA) 
{
  this=fcs
  
  if (is.na(staintable)) staintable=this$current.staintable
  
  if (is.na(index)) {
    printf("Warning: No file index chosen. Set to 1.")
    file=this$current.filenames[1]
    index=this$current.filetable[which(this$current.filetable[,2]==file),1]
  }
  
  vars=staintable[which(staintable[,1]==index),4]
  
  this$selected.vars = vars
  
  vars
}

##### Get data table from database
##### @param
# table         project name
# fileidx       file index
# stain         stain vector (optional)
fcs$getData <- function (
  table,
  fileidx,
  stain=NA,
  cofactor=NA) 
{
  this=fcs
  
  # get table from database
  data=dbGetQuery(this$conn,paste0("SELECT * FROM ",table," WHERE file_ID == '",fileidx,"'"))
  
  # cut column "file_ID" and ignore columns with NAs
  data$file_ID=NULL
  col.NA=NULL
  for (i in 1:ncol(data)) {
    if (any(is.na(data[,i]))) col.NA=c(col.NA,i)
  }
  if (!is.null(col.NA)) data=data[,-col.NA]
  
  this$current.project = table
  this$current.filenum = fileidx
  this$current.vars = this$selected.vars
  this$current.cutoffs = this$selected.cutoffs
  this$current.checks = this$selected.checks
  
  # set asinh cofactor to 1 if not set
  if ( is.na(cofactor) ) cofactor = 1
  
  data = asinh(data/cofactor) 
  
  if ( !is.na(stain) ) { 
    colnames(data) = stain 
  } else { 
    colnames(data) = this$selected.vars 
  }
  
  this$data = data
  
  printf("w: do getData: %s cells from table='%s' with fileidx=%s and cofactor=%s",nrow(data),table,fileidx,cofactor)
  
  data
}

##### Get cutoffs from current file
##### @param
# staintable    stain vector (optional)
# index         file index (optional)   
fcs$getCutoffs <- function (
  staintable=NA,
  index=NA) 
{
  this=fcs
  
  if (is.na(staintable)) staintable=this$current.staintable
  
  if (is.na(index)) {
    printf("Warning: No file index chosen. Set to 1.")
    file=this$current.filenames[1]
    index=this$current.filetable[which(this$current.filetable[,2]==file),1]
  }
  
  cutoffs=staintable[which(staintable[,1]==index),5]
  
  this$selected.cutoffs = cutoffs
  
  cutoffs
}

########################### Data import end


########################### Data visualization
##### add cutoff lines
##### @param
# cutoffs       cutoff vector (x,y)
fcs$addProdline <- function (cutoffs=c(0,0)) {
  this=fcs
  # add production cutoff line if provided
  if ( cutoffs[1] > 0 ) {
    abline(v=cutoffs[1],col="darkgrey")
  }
  if ( cutoffs[2] > 0 ) {
    abline(h=cutoffs[2],col="darkgrey")
  }
}


########################### Data manipulation
##### Calculate triplot matrices as a list
##### @param
# temp.data   tempory data frame with 3 columns v1,v2,v3
fcs$calc_triplot_mat <- function (
  temp.data,
  calc.meth,
  meta = NA,
  xy.range = NA,
  size.bin = 0.2,
  min.cells = 5) 
{
  is.even <- function(x) x %% 2 == 0
  
  brackets.open = c("(","[")
  
  ### cell count
  total.cells = nrow(temp.data)
  
  ### min/max values for x and y axis IF NOT SET
  if ( is.na(xy.range) ) {
    min.x = floor(min(temp.data[,1])*10)/10
    min.y = floor(min(temp.data[,2])*10)/10
    
    max.x = ceiling(max(temp.data[,1])*10)/10
    max.y = ceiling(max(temp.data[,2])*10)/10
  } else {
    min.x = min.y = xy.range[1]
    max.x = max.y = xy.range[2]
  }
  
  ### ranges need to be even!
  # only with bin size = 0.2
  if (size.bin == 0.2) {
    if (!is.even(min.x*10)) min.x = min.x - 0.1
    if (!is.even(min.y*10)) min.y = min.y - 0.1
    if (!is.even(max.x*10)) max.x = max.x +0.1
    if (!is.even(max.y*10)) max.y = max.y +0.1
  }
  
  ### initialize bin construct
  fX=cut(temp.data[,1],breaks=seq(min.x,max.x,by=size.bin),include.lowest=TRUE,dig.lab=5)
  fY=cut(temp.data[,2],breaks=seq(min.y,max.y,by=size.bin),include.lowest=TRUE,dig.lab=5)
  tab=table(fX,fY)
  seq.bin.y = seq(min.x,max.y-size.bin,by=size.bin)
  seq.bin.x = seq(min.y,max.x-size.bin,by=size.bin)
  colnames(tab) = seq.bin.y
  rownames(tab) = seq.bin.x
  ### combine X and Y position together
  fXY=as.factor(paste(fX,fY))
  
  ############ find minimum and maximum x where bins are displayed
  ### max.bin.y
  for ( max.bin.y in ncol(tab):1 ) {
    if ( any(tab[,max.bin.y]>=min.cells) ) break 
  }
  if ( max.bin.y<length(seq.bin.y) ) {
    max.bin.y = seq.bin.y[max.bin.y+1]
  } else {
    max.bin.y = max.y
  }
  ### min.bin.y
  for ( min.bin.y in 1:ncol(tab) ) {
    if ( any(tab[,min.bin.y]>=min.cells) ) break 
  }
  min.bin.y = seq.bin.y[min.bin.y]
  ############ find minimum and maximum y where bins are displayed
  ### max.bin.x
  for ( max.bin.x in nrow(tab):1 ) {
    if ( any(tab[max.bin.x,]>=min.cells) ) break 
  }
  if ( max.bin.x<length(seq.bin.x) ) {
    max.bin.x = seq.bin.x[max.bin.x+1]
  } else {
    max.bin.x = max.x
  }
  ### min.bin.x
  for ( min.bin.x in 1:nrow(tab) ) {
    if ( any(tab[min.bin.x,]>=min.cells) ) break 
  }
  min.bin.x = seq.bin.x[min.bin.x]
  
  ### group cells which calls into bin together and calculate bin method
  my.lengths = aggregate(temp.data[,3],by=list(fXY),length)
  
  if (calc.meth == "mean") {
    ### mean here is scaled mean (bin factor mean from 0-9)
    my.calc = aggregate(temp.data[,3], by=list(fXY), mean)
  } else if (calc.meth == "median") {
    my.calc = aggregate(temp.data[,3], by=list(fXY), median)
  } else if (calc.meth == "sd") {
    my.calc = aggregate(temp.data[,3], by=list(fXY), sd)
  } else if (calc.meth == "zrange") {
    my.calc = aggregate(temp.data[,3], by=list(fXY), function(z) {
      return(diff(range(max(z),min(z))))
    })
  }
  my.calc = cbind(my.calc, bin.cells = my.lengths$x)
  
  ### bin positions where we have at least minimum cells
  idx.min.cells = which(my.calc$bin.cells >= min.cells)
  ### min/max of range var colvec[v3]
  min.v3.range = floor(min(my.calc[idx.min.cells,'x'])*10)/10
  max.v3.range = ceiling(max(my.calc[idx.min.cells,'x'])*10)/10
  ### get bin steps of range var colvec[v3]
  step = round(diff(range(max.v3.range,min.v3.range))/10,2) 
  steps = seq(min.v3.range,max.v3.range, by=step)
  ### bin factor
  my.calc.fac = cut(my.calc$x, breaks=steps, labels=0:9, include.lowest=TRUE)
  
  ### combine length and bin factors
  my.calc = cbind(my.calc, fac=as.numeric(my.calc.fac)-1)
  
  mat = tab
  ### divide bin structure into quadrants
  for (x in rownames(tab)) {
    for (y in colnames(tab)) {
      ### only use bins with enough cells
      if ( tab[x,y] >= min.cells ) {
        ### create factor
        brackets.idx.x = brackets.idx.y = 1
        if (x==0) brackets.idx.x = 2
        if (y==0) brackets.idx.y = 2
        
        fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+size.bin,
                               '] ',brackets.open[brackets.idx.y],y,',',as.numeric(y)+size.bin,']',sep=''))
        
        idx.bin = which(as.character(fact)==as.character(my.calc$Group.1))
        fac.bin = my.calc[idx.bin,'fac']
        mat[x,y] = fac.bin
      } else {
        mat[x,y] = 'NA'
      }
    }
  }
  ranges=c(min.bin.x,max.bin.x,min.bin.y,max.bin.y,min.v3.range,max.v3.range)
  
  return(list(
    sample=meta$sample_id,
    condition=meta$condition,
    calc=calc.meth,
    markers=colnames(temp.data),
    "ranges(xmin,xmax,ymin,ymax,zmin,zmax)" = ranges,
    mat=mat
  )
  )
}

##### Calculate triplot matrices as a list in lappy
##### @param
# temp.data   tempory data frame with 3 columns v1,v2,v3
fcs$calc_triplot_mat_lapply <- function (
  temp.data,
  calc.meth,
  xy.range=NA) 
{
  is.even <- function(x) x %% 2 == 0
  
  global.min = 0
  global.max = 11.2
  xy.range = c(global.min,global.max)
  
  size.bin = 0.2
  min.cells = 5
  brackets.open = c("(","[")
  
  ### cell count
  total.cells = nrow(temp.data)
  
  ### min/max values for x and y axis IF NOT SET
  if ( is.na(xy.range) ) {
    min.x = floor(min(temp.data[,1])*10)/10
    min.y = floor(min(temp.data[,2])*10)/10
    
    max.x = ceiling(max(temp.data[,1])*10)/10
    max.y = ceiling(max(temp.data[,2])*10)/10
  } else {
    min.x = min.y = xy.range[1]
    max.x = max.y = xy.range[2]
  }
  
  ### ranges need to be even!
  # only with bin size = 0.2
  if (size.bin == 0.2) {
    if (!is.even(min.x*10)) min.x = min.x - 0.1
    if (!is.even(min.y*10)) min.y = min.y - 0.1
    if (!is.even(max.x*10)) max.x = max.x +0.1
    if (!is.even(max.y*10)) max.y = max.y +0.1
  }
  
  ### initialize bin construct
  fX=cut(temp.data[,1],breaks=seq(min.x,max.x,by=size.bin),include.lowest=TRUE,dig.lab=5)
  fY=cut(temp.data[,2],breaks=seq(min.y,max.y,by=size.bin),include.lowest=TRUE,dig.lab=5)
  tab=table(fX,fY)
  seq.bin.y = seq(min.x,max.y-size.bin,by=size.bin)
  seq.bin.x = seq(min.y,max.x-size.bin,by=size.bin)
  colnames(tab) = seq.bin.y
  rownames(tab) = seq.bin.x
  ### combine X and Y position together
  fXY=as.factor(paste(fX,fY))
  
  ############ find minimum and maximum x where bins are displayed
  ### max.bin.y
  for ( max.bin.y in ncol(tab):1 ) {
    if ( any(tab[,max.bin.y]>=min.cells) ) break 
  }
  if ( max.bin.y<length(seq.bin.y) ) {
    max.bin.y = seq.bin.y[max.bin.y+1]
  } else {
    max.bin.y = max.y
  }
  ### min.bin.y
  for ( min.bin.y in 1:ncol(tab) ) {
    if ( any(tab[,min.bin.y]>=min.cells) ) break 
  }
  min.bin.y = seq.bin.y[min.bin.y]
  ############ find minimum and maximum y where bins are displayed
  ### max.bin.x
  for ( max.bin.x in nrow(tab):1 ) {
    if ( any(tab[max.bin.x,]>=min.cells) ) break 
  }
  if ( max.bin.x<length(seq.bin.x) ) {
    max.bin.x = seq.bin.x[max.bin.x+1]
  } else {
    max.bin.x = max.x
  }
  ### min.bin.x
  for ( min.bin.x in 1:nrow(tab) ) {
    if ( any(tab[min.bin.x,]>=min.cells) ) break 
  }
  min.bin.x = seq.bin.x[min.bin.x]
  
  ### group cells which calls into bin together and calculate bin method
  my.lengths = aggregate(temp.data[,3],by=list(fXY),length)
  
  if (calc.meth == "mean") {
    ### mean here is scaled mean (bin factor mean from 0-9)
    my.calc = aggregate(temp.data[,3], by=list(fXY), mean)
  } else if (calc.meth == "median") {
    my.calc = aggregate(temp.data[,3], by=list(fXY), median)
  } else if (calc.meth == "sd") {
    my.calc = aggregate(temp.data[,3], by=list(fXY), sd)
  } else if (calc.meth == "zrange") {
    my.calc = aggregate(temp.data[,3], by=list(fXY), function(z) {
      return(diff(range(max(z),min(z))))
    })
  }
  my.calc = cbind(my.calc, bin.cells = my.lengths$x)
  
  ### bin positions where we have at least minimum cells
  idx.min.cells = which(my.calc$bin.cells >= min.cells)
  ### min/max of range var colvec[v3]
  min.v3.range = floor(min(my.calc[idx.min.cells,'x'])*10)/10
  max.v3.range = ceiling(max(my.calc[idx.min.cells,'x'])*10)/10
  ### get bin steps of range var colvec[v3]
  step = round(diff(range(max.v3.range,min.v3.range))/10,2) 
  steps = seq(min.v3.range,max.v3.range, by=step)
  ### bin factor
  my.calc.fac = cut(my.calc$x, breaks=steps, labels=0:9, include.lowest=TRUE)
  
  ### combine length and bin factors
  my.calc = cbind(my.calc, fac=as.numeric(my.calc.fac)-1)
  
  mat = tab
  ### divide bin structure into quadrants
  for (x in rownames(tab)) {
    for (y in colnames(tab)) {
      ### only use bins with enough cells
      if ( tab[x,y] >= min.cells ) {
        ### create factor
        brackets.idx.x = brackets.idx.y = 1
        if (x==0) brackets.idx.x = 2
        if (y==0) brackets.idx.y = 2
        
        fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+size.bin,
                               '] ',brackets.open[brackets.idx.y],y,',',as.numeric(y)+size.bin,']',sep=''))
        
        idx.bin = which(as.character(fact)==as.character(my.calc$Group.1))
        fac.bin = my.calc[idx.bin,'fac']
        mat[x,y] = fac.bin
      } else {
        mat[x,y] = 'NA'
      }
    }
  }
  
  return(mat)
}


##### NEW NEW NEW Calculate triplot quadrants function
##### @param
# temp.data     temporary data with 3 columns v1,v2,v3
# calc.meth     calculation method to use [mean,sd,median,zrange,range.quad,range.quad_scaled,var_quad]
# prod.cutoff   cutoff to divide neg/pos cells, here to divide quadrants. Be aware! Changes with different cofactors for asinh transformation
# size.bin      bin size
# min.cells     minimum number of cells per bin to consider this a countable bin
fcs$calc_triplot_quadrant <- function (
  temp.data, 
  calc.meth,
  prod.cutoff = NA,
  size.bin = 0.2, 
  min.cells = 5,
  min.bin.count = NA) 
{
  brackets.open = c("(","[")
    
  is.even <- function(x) x %% 2 == 0
  
  #min.x = min.y = 0
  ### min value for x and y axis
  min.x = floor(min(temp.data[,1])*10)/10
  min.y = floor(min(temp.data[,2])*10)/10
  
  ### max value for x and y axis
  max.x = ceiling(max(temp.data[,1])*10)/10
  max.y = ceiling(max(temp.data[,2])*10)/10
  
  ### ranges need to be even!
  # only with bin size = 0.2
  if (size.bin == 0.2) {
    if (!is.even(min.x*10)) min.x = min.x - 0.1
    if (!is.even(min.y*10)) min.y = min.y - 0.1
    if (!is.even(max.x*10)) max.x = max.x +0.1
    if (!is.even(max.y*10)) max.y = max.y +0.1
  }
  
  ### initialize bin construct
  fX=cut(temp.data[,1],breaks=seq(min.x,max.x,by=size.bin),include.lowest=TRUE,dig.lab=5)
  fY=cut(temp.data[,2],breaks=seq(min.y,max.y,by=size.bin),include.lowest=TRUE,dig.lab=5)
  tab=table(fX,fY)
  seq.bin.y = seq(min.y,max.y-size.bin,by=size.bin)
  seq.bin.x = seq(min.x,max.x-size.bin,by=size.bin)
  
  ### if minimum bin count is set
  if (!is.na(min.bin.count)) {
    
    ### if triploT has not enough bins, quadrant values will be NA, otherwise continue
    if (length(which(tab>min.cells)) < min.bin.count) {
      vec.q1 = vec.q2 = vec.q3 = vec.q4 = NA
    } else {
      ############ find minimum and maximum x where bins are displayed
      ### max.bin.y
      for ( max.bin.y in ncol(tab):1 ) {
        if ( any(tab[,max.bin.y]>=min.cells) ) break 
      }
      if ( max.bin.y<length(seq.bin.y) ) {
        max.bin.y = seq.bin.y[max.bin.y+1]
      } else {
        max.bin.y = max.y
      }
      ### min.bin.y
      for ( min.bin.y in 1:ncol(tab) ) {
        if ( any(tab[,min.bin.y]>=min.cells) ) break 
      }
      min.bin.y = seq.bin.y[min.bin.y]
      ############ find minimum and maximum y where bins are displayed
      ### max.bin.x
      for ( max.bin.x in nrow(tab):1 ) {
        if ( any(tab[max.bin.x,]>=min.cells) ) break 
      }
      if ( max.bin.x<length(seq.bin.x) ) {
        max.bin.x = seq.bin.x[max.bin.x+1]
      } else {
        max.bin.x = max.x
      }
      ### min.bin.x
      for ( min.bin.x in 1:nrow(tab) ) {
        if ( any(tab[min.bin.x,]>=min.cells) ) break 
      }
      min.bin.x = seq.bin.x[min.bin.x]
      
      ########### mean values for x and y where bins are displayed to divide into 4 quadrants
      # if cutoff is not set
      if (is.na(prod.cutoff)) {
        mean.bin.x = min.bin.x + (max.bin.x-min.bin.x)/2
        mean.bin.y = min.bin.y + (max.bin.y-min.bin.y)/2
      } else {
        mean.bin.x = mean.bin.y = prod.cutoff
      }
      
      ### if there is not enough cells
      if ( (min.bin.x > max.bin.x) | (min.bin.y > max.bin.y) ) {
        vec.q1 = vec.q2 = vec.q3 = vec.q4 = NA
      } else {
        ### make NEW bin construct with NEW minimum and maximum x/y where bins are displayed
        fX=cut(temp.data[,1],breaks=seq(min.bin.x,max.bin.x,by=size.bin),include.lowest=TRUE,dig.lab=5)
        fY=cut(temp.data[,2],breaks=seq(min.bin.y,max.bin.y,by=size.bin),include.lowest=TRUE,dig.lab=5)
        tab=table(fX,fY)
        
        if ( (max.bin.y - size.bin) > min.bin.y ) colnames(tab) = seq(min.bin.y,max.bin.y-size.bin,by=size.bin)
        else colnames(tab) = min.bin.y
        
        if ( (max.bin.x - size.bin) > min.bin.x ) rownames(tab) = seq(min.bin.x,max.bin.x-size.bin,by=size.bin)
        else rownames(tab) = min.bin.x
        # seq.bin.y = seq(min.bin.y,max.bin.y-size.bin,by=size.bin)
        # seq.bin.x = seq(min.bin.x,max.bin.x-size.bin,by=size.bin)
        # colnames(tab) = seq.bin.y
        # rownames(tab) = seq.bin.x
        ### combine X and Y position together
        fXY=as.factor(paste(fX,fY))
        
        ### group cells which calls into bin together and calculate bin method
        my.lengths = aggregate(temp.data[,3],by=list(fXY),length)
        
        if (calc.meth == "median") {
          my.calc = aggregate(temp.data[,3], by=list(fXY), median)
        } else if (calc.meth == "sd") {
          my.calc = aggregate(temp.data[,3], by=list(fXY), sd)
        } else if (calc.meth == "zrange") {
          my.calc = aggregate(temp.data[,3], by=list(fXY), function(z) {
            return(diff(range(max(z),min(z))))
          })
        } else if (calc.meth == "variance") {
          my.calc = aggregate(temp.data[,3], by=list(fXY), var)
        } else {   
          # if (calc.meth == "mean" | calc.meth == "range.quad" | calc.meth == "range.quad_scaled") 
          my.calc = aggregate(temp.data[,3], by=list(fXY), mean)
        } 
        
        my.calc = cbind(my.calc, bin.cells = my.lengths$x)
        
        ### bin positions where we have at least minimum cells
        idx.min.cells = which(my.calc$bin.cells >= min.cells)
        ### min/max of range var colvec[v3]
        min.v3.range = floor(min(my.calc[idx.min.cells,'x'])*10)/10
        max.v3.range = ceiling(max(my.calc[idx.min.cells,'x'])*10)/10
        ### get bin steps of range var colvec[v3]
        step = round(diff(range(max.v3.range,min.v3.range))/10,2) 
        steps = seq(min.v3.range,max.v3.range, by=step)
        ### bin factor
        ### e.g. mean here is scaled mean (bin factor mean from 0-9)
        my.calc.fac = cut(my.calc$x, breaks=steps, labels=0:9, include.lowest=TRUE)
        
        ### combine length and bin factors
        my.calc = cbind(my.calc, fac=as.numeric(my.calc.fac)-1)
        
        ### initialize quadrant vectors
        #vec.q1 = 
        vec.q2 = vec.q3 = vec.q4 = vector()
        
        ### divide bin structure into quadrants
        for (x in rownames(tab)) {
          for (y in colnames(tab)) {
            ### only use bins with enough cells
            if ( tab[x,y] >= min.cells ) {
              ### create factor
              brackets.idx.x = brackets.idx.y = 1
              if (x==0) brackets.idx.x = 2
              if (y==0) brackets.idx.y = 2
              
              fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+size.bin,
                                     '] ',brackets.open[brackets.idx.y],y,',',as.numeric(y)+size.bin,']',sep=''))
              
              idx.bin = which(as.character(fact)==as.character(my.calc$Group.1))
              
              if (calc.meth == "range.quad") {
                fac.bin = my.calc[idx.bin,'x']
              } else {
                fac.bin = my.calc[idx.bin,'fac']
              }
              
              #if ( (as.numeric(x) < mean.bin.x) & (as.numeric(y) < mean.bin.y) ) {
              ### Q1
              #  vec.q1 = c(vec.q1, fac.bin)
              #} else 
              if ( (as.numeric(x) < mean.bin.x) & (as.numeric(y) >= mean.bin.y) ) {
                ### Q2
                vec.q2 = c(vec.q2, fac.bin)
              } else if ( (as.numeric(x) >= mean.bin.x) & (as.numeric(y) >= mean.bin.y) ) {
                ### Q3
                vec.q3 = c(vec.q3, fac.bin)
              } else if ( (as.numeric(x) >= mean.bin.x) & (as.numeric(y) < mean.bin.y) ) {
                ### Q4
                vec.q4 = c(vec.q4, fac.bin)
              }
            }
          }
        }
      }
    }
  }
  
  ### return range of every quadrants listed bins with minimum cells and used calculation method
  if (calc.meth == "range.quad") {
    return( c(#round(diff(range(vec.q1)),2),
      round(diff(range(vec.q2)),2),
      round(diff(range(vec.q3)),2),
      round(diff(range(vec.q4)),2)
    ))
  ### return variance of every quadrants listed bins with minimum cells and used calculation method
  } else if (calc.meth == "var_quad") {
    return( c(#round(mean(vec.q1),2),
      round(var(vec.q2),2),
      round(var(vec.q3),2),
      round(var(vec.q4),2) 
    ))
  ### return mean of every quadrants listed bins with minimum cells and used calculation method
  } else {
    return( c(round(mean(vec.q1),2),
      round(mean(vec.q2),2),
      round(mean(vec.q3),2),
      round(mean(vec.q4),2) 
    ))
  }
      
}


##### Calculate triplot quadrants function by MAX-POOLING
##### @param
# temp.data     temporary data with 3 columns v1,v2,v3 [data.frame]
# calc.meth     calculation method to use ["mean","sd","median","zrange"]
# prod.cutoff   cutoff to divide neg/pos cells, here to divide quadrants. Be aware! Changes with different cofactors for asinh transformation [var]
# size.bin      bin size [var]
# min.cells     minimum number of cells per bin to consider this a countable bin [var]
# min.bin.count minimum number of bins in this 3 combinatoric constellation [var]
# minimize.area minimize bin area [TRUE/FALSE]
fcs$calc_triplot_maxpool <- function (
  temp.data, 
  calc.meth,
  prod.cutoff = NA,
  size.bin = 0.2, 
  min.cells = 5,
  maxp.size = 3,
  stride = 3,
  min.bin.count = 50,
  minimize.area = FALSE,
  xy.range = c(0,9,0,9)) 
{
  brackets.open = c("(","[")
  
  is.even <- function(x) x %% 2 == 0
  
  
  ### min/max values for x and y axis IF NOT SET
  if ( minimize.area ) {
    min.x = floor(min(temp.data[,1])*10)/10
    min.y = floor(min(temp.data[,2])*10)/10
    
    max.x = ceiling(max(temp.data[,1])*10)/10
    max.y = ceiling(max(temp.data[,2])*10)/10
  } else {
    min.x = min.y = xy.range[1]
    max.x = max.y = xy.range[2]
  }

  ### ranges need to be even!
  # only with bin size = 0.2
  if (size.bin == 0.2) {
    if (!is.even(min.x*10)) min.x = min.x - 0.1
    if (!is.even(min.y*10)) min.y = min.y - 0.1
    if (!is.even(max.x*10)) max.x = max.x +0.1
    if (!is.even(max.y*10)) max.y = max.y +0.1
  }
  
  ### initialize bin construct
  fX=cut(temp.data[,1],breaks=seq(min.x,max.x,by=size.bin),include.lowest=TRUE,dig.lab=5)
  fY=cut(temp.data[,2],breaks=seq(min.y,max.y,by=size.bin),include.lowest=TRUE,dig.lab=5)
  tab=table(fX,fY)
  seq.bin.y = seq(min.y,max.y-size.bin,by=size.bin)
  seq.bin.x = seq(min.x,max.x-size.bin,by=size.bin)
  
  
  ### if triploT has not enough bins, quadrant values will be NA, otherwise continue
  if ( length(which(tab>min.cells)) < min.bin.count ) {
    q1 = q2 = q3 = q4 = NA
  } else {
    
    
    ### do not include this, try without minimizing bin area
    if ( minimize.area ) {
      ############ find minimum and maximum x where bins are displayed
      ### max.bin.y
      for ( max.bin.y in ncol(tab):1 ) {
        if ( any(tab[,max.bin.y]>=min.cells) ) break 
      }
      if ( max.bin.y<length(seq.bin.y) ) {
        max.bin.y = seq.bin.y[max.bin.y+1]
      } else {
        max.bin.y = max.y
      }
      ### min.bin.y
      for ( min.bin.y in 1:ncol(tab) ) {
        if ( any(tab[,min.bin.y]>=min.cells) ) break 
      }
      min.bin.y = seq.bin.y[min.bin.y]
      ############ find minimum and maximum y where bins are displayed
      ### max.bin.x
      for ( max.bin.x in nrow(tab):1 ) {
        if ( any(tab[max.bin.x,]>=min.cells) ) break 
      }
      if ( max.bin.x<length(seq.bin.x) ) {
        max.bin.x = seq.bin.x[max.bin.x+1]
      } else {
        max.bin.x = max.x
      }
      ### min.bin.x
      for ( min.bin.x in 1:nrow(tab) ) {
        if ( any(tab[min.bin.x,]>=min.cells) ) break 
      }
      min.bin.x = seq.bin.x[min.bin.x]
      
      ########### mean values for x and y where bins are displayed to divide into 4 quadrants
      # if cutoff is not set
      if (is.na(prod.cutoff)) {
        mean.bin.x = min.bin.x + (max.bin.x-min.bin.x)/2
        mean.bin.y = min.bin.y + (max.bin.y-min.bin.y)/2
      } else {
        mean.bin.x = mean.bin.y = prod.cutoff
      }
      
      ### make NEW bin construct with NEW minimum and maximum x/y where bins are displayed
      fX=cut(temp.data[,1],breaks=seq(min.bin.x,max.bin.x,by=size.bin),include.lowest=TRUE,dig.lab=5)
      fY=cut(temp.data[,2],breaks=seq(min.bin.y,max.bin.y,by=size.bin),include.lowest=TRUE,dig.lab=5)
      tab=table(fX,fY)
      
      if ( (max.bin.y - size.bin) > min.bin.y ) {
        colnames(tab) = seq(min.bin.y,max.bin.y-size.bin,by=size.bin)
      } else {
        colnames(tab) = min.bin.y
      }
      
      if ( (max.bin.x - size.bin) > min.bin.x ) {
        rownames(tab) = seq(min.bin.x,max.bin.x-size.bin,by=size.bin)
      } else {
        rownames(tab) = min.bin.x
      }
    }
    
    
    # seq.bin.y = seq(min.bin.y,max.bin.y-size.bin,by=size.bin)
    # seq.bin.x = seq(min.bin.x,max.bin.x-size.bin,by=size.bin)
    # colnames(tab) = seq.bin.y
    # rownames(tab) = seq.bin.x
    ### combine X and Y position together
    fXY=as.factor(paste(fX,fY))
    
    ### group cells which calls into bin together and calculate bin method
    my.lengths = aggregate(temp.data[,3],by=list(fXY),length)
    
    if (calc.meth == "median") {
      my.calc = aggregate(temp.data[,3], by=list(fXY), median)
    } else if (calc.meth == "sd") {
      my.calc = aggregate(temp.data[,3], by=list(fXY), sd)
    } else if (calc.meth == "zrange") {
      my.calc = aggregate(temp.data[,3], by=list(fXY), function(z) {
        return(diff(range(max(z),min(z))))
      })
    } else if (calc.meth == "variance") {
      my.calc = aggregate(temp.data[,3], by=list(fXY), var)
    } else {   
      # if (calc.meth == "mean" | calc.meth == "range.quad" | calc.meth == "range.quad_scaled") 
      my.calc = aggregate(temp.data[,3], by=list(fXY), mean)
    }
    my.calc$x = round(my.calc$x,2)
    my.calc = cbind(my.calc, bin.cells = my.lengths$x)
    
    ### bin positions where we have at least minimum cells
    idx.min.cells = which(my.calc$bin.cells >= min.cells)
    ### min/max of range var colvec[v3]
    min.v3.range = floor(min(my.calc[idx.min.cells,'x'])*10)/10
    max.v3.range = ceiling(max(my.calc[idx.min.cells,'x'])*10)/10
    ### get bin steps of range var colvec[v3]
    step = round(diff(range(max.v3.range,min.v3.range))/10,2) 
    steps = seq(min.v3.range,max.v3.range, by=step)
    ### bin factor
    ### e.g. mean here is scaled mean (bin factor mean from 0-9)
    my.calc.fac = cut(my.calc$x, breaks=steps, labels=0:9, include.lowest=TRUE)
    
    ### combine length and bin factors
    my.calc = cbind(my.calc, fac=as.numeric(my.calc.fac)-1)
    
    
    
    
    ### create calc matrix
    tab.calc = matrix( data = my.calc$x,
                       nrow = nrow(tab),
                       ncol = ncol(tab)
    )
    ### set NAs where minimum amount of cells are note reached
    tab.calc[which(tab<min.cells)] = NA
    rownames(tab.calc) = rownames(tab)
    colnames(tab.calc) = colnames(tab)
    
    
    ### decrease dimension
    it = 0
    while ( ncol(tab.calc)>(maxp.size+2*stride) & nrow(tab.calc)>(maxp.size+2*stride)) {
      it = it + 1
      printf("it=%s",it)
      
      tab.calc.row.stride = seq(1,nrow(tab.calc)-maxp.size,by=stride)
      tab.calc.col.stride = seq(1,ncol(tab.calc)-maxp.size,by=stride)
      xi = vector()
      for ( x in tab.calc.row.stride ) {
        print(x)
        yi = vector()
        for ( y in tab.calc.col.stride ) {
          printf("y=%s",y)
          print(tab.calc[x:(x+maxp.size-1),y:(y+maxp.size-1)])
          
          tmp = tab.calc[x:(x+maxp.size-1),y:(y+maxp.size-1)]
          if (any(is.na(tmp))) tmp = tmp[-which(is.na(tmp))]
          
          if (length(tmp)!=0) {
            yi = c(yi,max(tmp))
            printf("max=%s",max(tmp))
          } else {
            yi = c(yi,NA)
            print("NA")
          }
          
          
        }
        print(yi)
        xi=rbind(xi,yi)
        
        
      }
      tab.calc=xi
      
      print(tab.calc)
    }
    
    ### last max pooling step
    ## divide bin structure into 4 pools
    ## row == col in "real", acutally
    col.half = ceiling((ncol(tab.calc)/2))
    row.half = ceiling((nrow(tab.calc)/2))
    
    ### q1
    q1 = tab.calc[1:row.half,1:col.half]
    if (any(is.na(q1))) q1 = q1[-which(is.na(q1))]
    q1 = max(q1)
    ### q2
    q2 = tab.calc[(row.half+1):nrow(tab.calc),1:col.half]
    if (any(is.na(q2))) q2 = q2[-which(is.na(q2))]
    q2 = max(q2)                                 
    ### q3
    q3 = tab.calc[(row.half+1):nrow(tab.calc),(col.half+1):ncol(tab.calc)]
    if (any(is.na(q3))) q3 = q3[-which(is.na(q3))]
    q3 = max(q3)
    ### q4
    q4 = tab.calc[1:row.half,(col.half+1):ncol(tab.calc)]
    if (any(is.na(q4))) q4 = q4[-which(is.na(q4))]
    q4 = max(q4)
  }
  
  
  return(c(q1,q2,q3,q4))
  
}

##### OLD OLD OLD Calculate triplot quadrants function OLD OLD OLD
##### @param
# temp.data        temporary data frame with 3 parameters v1,v2,v3
# cutoffs     cutoff vector with 3 variables, default=(0,0,0)
# binSize     (asinh) size of bin, default=0.2
# mincells    minimum number of cells where bin would be colored, default=10
# checkCALC   calculation method, default="mean"
# data.origin original temporary data frame with 3 parameters, only applicable if checkCALC="mean+"
fcs$calc_triplot_quadrants <- function(
  temp.data,
  cutoffs=c(0,0,0),
  binSize=0.2,
  mincells=10,
  checkCALC="mean",
  data.origin=NA) 
{
  this=fcs
  
  
  
  prodcells.color="red"
  prodpluscells.color="chartreuse4"
  
  file=tclvalue(tkget(this$tkchoosefile))
  displayfile = this$shortenFilename(file)
  
  metadatafile = sprintf("%s_PRIvis_metadata.csv",this$current.project)
  if(!file.exists(metadatafile)) {
    header = c("sample","feat.X","feat.Y", "feat.Z", "calc"
               ,"q1.total","q2.total","q2.total","q2.total"
               ,"q1.prodcells","q2.prodcells","q3.prodcells","q4.prodcells"
               ,"q1.prodcellsplus","q2.prodcellsplus","q3.prodcellsplus","q4.prodcellsplus")
    write.table(t(header), metadatafile, sep = ",", row.names=F, col.names=F)
  }
  
  data = as.matrix(data)
  if (!is.na(data.origin)) data.origin = as.matrix(data.origin)
  
  # axes range
  xmin.val=as.numeric(tkget(this$minvalX))
  xmax.val=as.numeric(tkget(this$maxvalX))
  ymin.val=as.numeric(tkget(this$minvalY))
  ymax.val=as.numeric(tkget(this$maxvalY))
  min.mFI=as.double(tclvalue(this$vminMFI))
  max.mFI=as.double(tclvalue(this$vmaxMFI))
  
  # checkbutton options
  checkDYNRANGE = tclvalue(this$cbtdynRange)
  checkTRANS = tclvalue(this$rbtrans)
  if (checkTRANS =="") checkTRANS = tclvalue(this$rbtrans) = "asinh"
  checkCALC = tclvalue(this$rbcalc)
  if (checkCALC == "") checkCALC = tclvalue(this$rbcalc) = "MFI"
  if (checkCALC == "density") density=TRUE
  checkGRID = tclvalue(this$cbtshowGrid)
  
  if (this$working) printf("w: do bintriplot density=%s checkCALC=%s",density,checkCALC)
  if (checkTRANS=="asinh") options(scipen=-1)
  else options(scipen=999)
  
  # legend from blue to red
  if (bg) cols=rep("gray",12)
  else cols=this$col.rainbow
  
  # legend title
  tmp = unlist(strsplit(colnames(data)[3],"\\."))
  if ( length(tmp) > 1 ) {
    if ( cutoffs[3]>0 ) legend.title = sprintf("%s(%s)",tmp[1],cutoffs[3]) 
    else legend.title = tmp[1]
  } else {
    if ( cutoffs[3]>0 ) legend.title = sprintf("%s(%s)",colnames(data)[3],cutoffs[3]) 
    else legend.title = colnames(data)[3]
  }
  
  # boolean for only grey plot 
  # if MFI(+) mode and there are no colorful bins to display
  grey.label = TRUE
  my.LENGTHS = FALSE
  
  # set negative values of z-axis (colnum=3) to zero
  if (!density) data[which(data[,3]<0),3] = 0
  
  ncells=nrow(data)
  
  if ( ncells > 0 ) {
    
    ### tdata = cut cells which lie in plot area
    #tdata2 = data[-c(which(data[,1]<xmin.val),which(data[,2]<xmin.val)), ]
    tdata = data[which(data[,1]>=xmin.val & data[,2]>=ymin.val), ]
    
    ### construct bins
    fX=cut(tdata[,1],breaks=seq(xmin.val,xmax.val,by=binSize),include.lowest=TRUE,dig.lab=5)
    fY=cut(tdata[,2],breaks=seq(ymin.val,ymax.val,by=binSize),include.lowest=TRUE,dig.lab=5)
    tab=table(fX,fY)
    
    colnames(tab)=seq(ymin.val,ymax.val-binSize,by=binSize)
    rownames(tab)=seq(xmin.val,xmax.val-binSize,by=binSize)
    fXY=as.factor(paste(fX,fY))
    
    ### construct bins for MFI(+)
    if (checkCALC == "MFI(+)") {
      fX.origin=cut(data.origin[,1],breaks=seq(xmin.val,xmax.val,by=binSize),include.lowest=TRUE,dig.lab=5)
      fY.origin=cut(data.origin[,2],breaks=seq(ymin.val,ymax.val,by=binSize),include.lowest=TRUE,dig.lab=5)
      tab.origin=table(fX.origin,fY.origin)
      colnames(tab.origin)=seq(ymin.val,ymax.val-binSize,by=binSize)
      rownames(tab.origin)=seq(xmin.val,xmax.val-binSize,by=binSize)
    } else {
      tab.origin = tab
    }
    
    if (density) {
      # number of cells in bin
      my.calc=aggregate(tdata[,3],by=list(fXY),length)
    } else {
      # get means/median/freq
      if ( grepl("MFI",checkCALC) ) {
        my.calc=aggregate(tdata[,3],by=list(fXY),mean)
        #else if ( checkCALC == "medianFI" )  my.calc=aggregate(tdata[,3],by=list(fXY),median)
      } else if ( checkCALC == "SD" ) {
        my.calc=aggregate(tdata[,3],by=list(fXY),sd)
        cols=this$col.blackwhite
      } else if ( checkCALC == "SEM" ) {
        my.calc=aggregate(tdata[,3],by=list(fXY),function(x) {
          SEM = sd(x)/sqrt(length(x))
          # if normally distributed, 95,4 % of the cells should lie inside the interval mean +/- SEM
          #interval_min = mean(x) - SEM
          #interval_max = mean(x) + SEM 
          SEM
        })
        cols=this$col.blackwhite
      } else if ( checkCALC == "RSEM" ) {
        my.calc=aggregate(tdata[,3],by=list(fXY),function(x) {
          RSEM = sd(x)/sqrt(length(x))
          RSEM/mean(x)*100
        })
      } else if ( checkCALC == "freq" ) {
        my.calc = aggregate(tdata[,3],by=list(fXY),function(x) {
          y= round( 100 * length(which(x >= cutoffs[3])) / length(x))
          return(y)
        })
      } 
    }
    ### http://www.allgemeinmedizin.med.uni-goettingen.de/de/media/2008_Koschack_Standardabweichung_Standardfehler.pdf
    ##  Standardabweichung (SD) 
    # – ist eine Aussage über die Streuung der erhobenen Werte in einer Stichprobe
    # – ist nur wenig durch die Grösse der Stichprobe beeinflussbar
    # – hängt von der biologischen Variabilität ab
    # – ist ein beschreibendes Mass
    ## Standardfehler (SEM)
    # – ist eine Aussage über die „Genauigkeit“ des Mittelwerts in einer Stichprobe
    # – hängt von der Messgenauigkeit ab
    # – ist ein statistisches Mass
    # – steht in direktem Verhältnis zur Grösse der Stichprobe
    
    ### https://www.graphpad.com/guides/prism/6/statistics/index.htm?stat_semandsdnotsame.htm
    # It is easy to be confused about the difference between the standard deviation (SD) and the standard error of the mean (SEM). Here are the key differences:
    # •   The SD quantifies scatter — how much the values vary from one another.
    # •   The SEM quantifies how precisely you know the true mean of the population. It takes into account both the value of the SD and the sample size.
    # •   Both SD and SEM are in the same units -- the units of the data.
    # •   The SEM, by definition, is always smaller than the SD.
    # •   The SEM gets smaller as your samples get larger. This makes sense, because the mean of a large sample is likely to be closer to the true population mean than is the mean of a small sample. With a huge sample, you'll know the value of the mean with a lot of precision even if the data are very scattered.
    # •   The SD does not change predictably as you acquire more data. The SD you compute from a sample is the best possible estimate of the SD of the overall population. As you collect more data, you'll assess the SD of the population with more precision. But you can't predict whether the SD from a larger sample will be bigger or smaller than the SD from a small sample. (This is not strictly true. It is the variance -- the SD squared -- that doesn't change predictably, but the change in SD is trivial and much much smaller than the change in the SEM.)
    # Note that standard errors can be computed for almost any parameter you compute from data, not just the mean. The phrase "the standard error" is a bit ambiguous. The points above refer only to the standard error of the mean.
    
    
    my.lengths=aggregate(tdata[,3],by=list(fXY),length)
    
    my.LENGTHS = any(my.lengths[,2] >= mincells)
  }
  
  ### if there are bins to display
  if ( my.LENGTHS ) {
    # there are bins to plot, so set grey.label to FALSE
    grey.label = FALSE
    
    my.calc=cbind(my.calc,ncells=my.lengths$x)
    
    if ( checkCALC == "MEAN_SEM" & !density) {
      this$my.test = aggregate(tdata[,3],by=list(fXY),function(x) {
        SEM = sd(x)/sqrt(length(x))
        
        if ( SEM >= 0.5 ) {
          return(1)
        } else { 
          return (0) 
        }
        # if normally distributed, 95,4 % of the cells should lie inside the interval mean +/- 2*SEM
      })
    }
    
    # delete rows with NA!
    # my.calc.NA = my.calc[-grep('NA',my.calc$Group.1),]
    
    ### get steps for legend and plot in mode 'freq'
    decim=1
    #range=""
    col.minmax="black"
    if (checkCALC == "freq" & !density) {
      #min.legend=ymin.val + 3*binSize
      #max.legend=max.range=(ymax.val-ymin.val)/2 #- (ymax.val-ymin.val)/10
      #max.legend = max.range = ymin.val + 3*binSize + diff(c(ymin.val,ymax.val))/2
      max.range = ymin.val + 3*binSize + diff(c(ymin.val,ymax.val))/2
      
      label.steps=seq(0,100,by=10)
      
      # bin color factor
      my.calc.fac=cut(my.calc$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
      levels(my.calc.fac)=c(0,levels(my.calc.fac),11,12)
      
      decim=0
    } else if (checkCALC == "RSEM" & !density) {
      #min.legend=ymin.val
      #max.legend=max.range=(ymax.val-ymin.val)/3
      min.legend = ymin.val + 3*binSize
      max.legend = max.range = ymin.val + 3*binSize + diff(c(ymin.val,ymax.val))/3
      
      step=round(diff(range(max.legend,min.legend))/6,1)
      steps=seq(min.legend,max.legend,by=step)
      label.steps=seq(25,55,by=5)
      
      label.steps[7]=""
      label.steps[6]=">=50"
      label.steps[2]="<=25"
      label.steps[1]=""
      
      # bin color factor
      my.calc.fac=cut(my.calc$x,breaks=seq(0,50,by=5),labels=1:10,include.lowest=TRUE)
      levels(my.calc.fac)=c(0,levels(my.calc.fac),11,12)
      this$my.calc.fac2 = my.calc.fac
      for ( i in 1:length(my.calc.fac) ){
        if ( !is.na(my.calc$x[i]) ) {
          if (my.calc$x[i]>=50 & my.calc$ncells[i]>=mincells) {
            my.calc.fac[i] = 11
          }
        }
      }
      cols = this$col.blackred
      decim=0
    } else {
      # get steps for legend and plot
      if (density) {
        idx=which(my.calc$ncells>=mincells)
        #idx=idx[grep("NA",my.calc[idx,'Group.1'])]
        min.range=floor(min(my.calc[idx,'x'])*10)/10
        max.range=max(tab)
        #max.range = max(my.lengths[,2])
        decim=0
        #printf("maxrange=%s",max.range)
        
        if ( max.range < 200 ) col.minmax="red"
      } else if ( checkDYNRANGE=="1" ) {
        idx=which(my.calc$ncells>=mincells)
        #idx=idx[grep("NA",my.calc[idx,'Group.1'])]
        min.range=floor(min(my.calc[idx,'x'])*10)/10
        max.range=ceiling(max(my.calc[idx,'x'])*10)/10
        
        # if dynamic range is too small
        if ( grepl("MFI",checkCALC) & diff(c(min.range,max.range)) <= 0.5 ) col.minmax="red"
      } else {
        min.range=min.mFI
        max.range=max.mFI
        #range=sprintf("Manual MFI range: %0.1f-%0.1f",min.range,max.range)
      }
      
      # get steps
      step=round(diff(range(max.range,min.range))/10,2) 
      steps=seq(min.range,max.range,by=step)
      label.steps=steps[1:11]
      if (density | checkTRANS=="biex") {
        if(max.range>500) label.steps = round(label.steps,-2)
        else label.steps = round(label.steps)
        if (density) label.steps[1] = min.range
      } else if ( checkDYNRANGE != "1" ) {
        label.steps[1]=sprintf("<=%s",label.steps[1])
        label.steps[11]=sprintf(">=%s",label.steps[11])
      } else {
        label.steps = round(label.steps,1)
      }
      
      # bin color factor
      my.calc.fac=cut(my.calc$x,breaks=steps,labels=2:11,include.lowest=TRUE)
      
      levels(my.calc.fac)=c(0,levels(my.calc.fac),12)
      # if x < min.range
      my.calc.fac[which(my.calc$x<steps[1] & my.calc$ncells>=mincells)]=0
      # if x > max.range
      my.calc.fac[which(my.calc$x>steps[11] & my.calc$ncells>=mincells)]=11
    }
    my.calc=cbind(my.calc,fac=as.numeric(my.calc.fac)+1)
    
    #this$my.calc.fac = my.calc.fac
    this$my.calc = my.calc
    
    this$bincount = 0
    this$maxcells = 0
    
    ##### plot bins
    for (x in rownames(tab)) {
      for (y in colnames(tab)) {
        if ( tab[x,y]>=mincells ) {
          fact=as.factor(paste('(',x,',',as.numeric(x)+binSize,'] ','(',y,',',as.numeric(y)+binSize,']',sep=''))
          idx=which(as.character(fact)==as.character(my.calc$Group.1))
          rect(x,y,as.numeric(x)+binSize,as.numeric(y)+binSize,col=cols[my.calc[idx,'fac']],border=NA)
          
          this$bincount = this$bincount + 1                
          
          if (tab[x,y]>this$maxcells) this$maxcells=tab[x,y]            
        } else if ( checkCALC == "MFI(+)" & tab.origin[x,y] >= mincells) {
          rect(x,y,as.numeric(x)+binSize,as.numeric(y)+binSize,col="gray",border=NA)
        }
      }
    }
  }
}

########################### Data manipulation end




