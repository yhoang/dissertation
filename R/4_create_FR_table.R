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

folder = sprintf("table/",tissue,checkCALC,alphalist)

file.name = sprintf("signif_quadrants%s_it%s.txt",quad.cap,iterations)
output.quads = sprintf("%s_%s_signif_quadrants%s_it%s_table.csv",tissue,checkCALC,quad.cap,iterations)
output.quad.vars = sprintf("%s_%s_signif_quadrant%s_it%s_vars_table.csv",tissue,checkCALC,quad.cap,iterations)

### get protein names    
protein.names = readRDS("meta/protein_names")
xyz.vec=c("_X","_Y","_Z")

signif.table = read.table(file=paste0(folder,file.name), na.strings = "NA", fill =TRUE,stringsAsFactors=F)
names(signif.table)[1:5] = c("tissue","calculation_meth","sample_size","prediction_RMSE","seed")

#quads = signif.table[which(signif.table$prediction_RMSE<0.05 & signif.table$sample_size=="quads"),6:ncol(signif.table)]
quads = signif.table[which(signif.table$sample_size=="quads"),6:ncol(signif.table)]
quads = as.vector(t(quads))

quads.filtered = quads[!is.na(quads)]
quads.filtered = quads.filtered[-which(quads.filtered=="")]
quads.filtered = sort(table(quads.filtered),decreasing = T)
quads.filtered = as.data.frame(quads.filtered)

quads.filtered.freq = round(sapply(quads.filtered$Freq, function(x) {x/sum(quads.filtered$Freq)}),3)


#quads.coeff = signif.table[which(signif.table$prediction_RMSE<0.05 & signif.table$sample_size=="coeff"),6:ncol(signif.table)]
quads.coeff = signif.table[which(signif.table$sample_size=="coeff"),6:ncol(signif.table)]
quads.coeff = as.vector(t(quads.coeff))
names(quads.coeff) = quads
quads.coeff = quads.coeff[-which(names(quads.coeff)=="")]

coeff.vec = vector()
for ( j in 1:nrow(quads.filtered) ) {
  quad.idx = which(names(quads.coeff)==quads.filtered[j,1])
  coeff.vec = c(coeff.vec,round(mean(as.numeric(quads.coeff[quad.idx])),4))
}

triplot.df = data.frame( matrix(NA,
  nrow = nrow(quads.filtered),
  ncol = 4,
  byrow = TRUE,
), stringsAsFactors = FALSE
)
for (it in 1:nrow(quads.filtered)) {
  var.it.list = unlist(strsplit(as.character(quads.filtered[it,1]),"[.]"))
  len.var = length(var.it.list)
  var.pos = vector()
  if ( len.var==6 ) {
    # no proteins with "." 
    triplot.df[it,] = var.it.list[c(1:3,len.var)]
  } else {
    dot.tmp = 0
    for (pos in 1:(length(var.it.list)-2)) {
      if ( (var.it.list[pos] %in% protein.names) ) {
        var.pos = c(var.pos,paste0(var.it.list[pos]))
      } else if (dot.tmp == 1) {
        var.tmp = paste0(var.it.list[pos-1],".",var.it.list[pos])
        var.pos = c(var.pos,var.tmp)
        dot.tmp = 0
      } else {
        dot.tmp = 1
      }
    }
    #triplot.df[it,] = c(gsub("_X|_Y|_Z","",var.pos),var.it.list[len.var])
    triplot.df[it,] = c(var.pos,var.it.list[len.var])
  }
}

quads.df = cbind(rownames(quads.filtered),quads.filtered,quads.filtered.freq,coeff.vec,triplot.df)
colnames(quads.df) = c("rank","quadrants","counts","% counts","coefficient","X","Y","Z","Q")

write.table(quads.df,sep="\t",col.names=T,row.names=F,quote=F,file=paste0(folder,output.quads))
printf("Written in %s.",paste0(folder,output.quads))

##########################
### count protein in protein combination of quadrants 
names(triplot.df) = c("X","Y","Z","Q")


###### single counts
## X
var.count.X = table(triplot.df$X)
var.count.X = as.data.frame(var.count.X)
rownames(var.count.X)=var.count.X[,1]
colnames(var.count.X) = c("marker_X","count_X")
var.count.X = var.count.X[order(var.count.X$count_X,decreasing = T),]
freq.X = round(var.count.X$count_X/sum(var.count.X$count_X),2)
var.count.X = cbind(var.count.X, freq_X = round(var.count.X$count_X/sum(var.count.X$count_X),3))
## Y
var.count.Y = table(triplot.df$Y)
var.count.Y = as.data.frame(var.count.Y)
rownames(var.count.Y)=var.count.Y[,1]
colnames(var.count.Y) = c("marker_Y","count_Y")
var.count.Y = var.count.Y[order(var.count.Y$count_Y,decreasing = T),]
var.count.Y = cbind(var.count.Y,freq_Y = round(var.count.Y$count_Y/sum(var.count.Y$count_Y),3))
## Z
var.count.Z = table(triplot.df$Z)
var.count.Z = as.data.frame(var.count.Z)
rownames(var.count.Z)=var.count.Z[,1]
colnames(var.count.Z) = c("marker_Z","count_Z")
var.count.Z = var.count.Z[order(var.count.Z$count_Z,decreasing = T),]
var.count.Z = cbind(var.count.Z,freq_Z = round(var.count.Z$count_Z/sum(var.count.Z$count_Z),3))

###### combined counts
## XY
var.count.XY = merge(var.count.X,var.count.Y,by="row.names",all=T)
rownames(var.count.XY) = var.count.XY$Row.names
var.count.XY = cbind(var.count.XY,marker_XY = var.count.XY$Row.names)
var.count.XY = var.count.XY[,-1]
var.count.XY[is.na(var.count.XY)] = 0
var.count.XY = cbind(var.count.XY, count_XY=var.count.XY$count_X+var.count.XY$count_Y)
var.count.XY = var.count.XY[order(var.count.XY$count_XY,decreasing = T),]
var.count.XY = cbind(var.count.XY,freq_XY = round(var.count.XY$count_XY/sum(var.count.XY$count_XY),3))
colnames(var.count.XY)
## XYZ
var.count.all = merge(var.count.XY,var.count.Z,by="row.names",all=T)
rownames(var.count.all) = var.count.all$marker_X = var.count.all$marker_Y = var.count.all$marker_Z = var.count.all$Row.names
var.count.all = cbind(var.count.all,marker_all = var.count.all$Row.names)
var.count.all = var.count.all[,-1]
var.count.all[is.na(var.count.all)] = 0
var.count.all = cbind(var.count.all, count_all=var.count.all$count_XY+var.count.all$count_Z)
var.count.all = var.count.all[order(var.count.all$count_all,decreasing = T),]
var.count.all = cbind(var.count.all, freq_all = round(var.count.all$count_all/sum(var.count.all$count_all),3))

var.count.all = cbind(var.count.all[order(var.count.all$count_X,decreasing = T),1:3],   # order to marker_X
                      var.count.all[order(var.count.all$count_Y,decreasing = T),4:6],   # order to marker_Y
                      var.count.all[order(var.count.all$count_XY,decreasing = T),7:9],   # order to marker_XY
                      var.count.all[order(var.count.all$count_Z,decreasing = T),10:12],   # order to marker_Z
                      var.count.all[order(var.count.all$count_all,decreasing = T),13:15])   # order to marker_Z

### write only counts of proteins from significant quadrants
# create file and header
write.table(var.count.all,sep="\t", col.names=T,row.names=F,quote=F,
            file=paste0(folder,output.quad.vars))

printf("Written in %s.",paste0(folder,output.quad.vars))
