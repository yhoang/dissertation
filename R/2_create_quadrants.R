#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
rm(list = ls())
source("R/spitzer.r")


meta = new.env()
checkCALC = "range.quad"


cl <- makeCluster(3)
cofactor = 0.1
min.cells = 5
bin.min.count = 400


## ----metadata, cache=TRUE, results="hide"--------------------------------
this=fcs
metadata = read.csv(file="meta/Spitzer_blood.csv")
metadata$condition =factor(metadata$condition, levels =  
                             c("Bl_CD1_d3","Bl_aPD1_d3","Bl_B6_d3","Bl_untr_d3"))

# dunkel
color_conditions = c("#006837","#b2182b","#2166ac","#8c510a")
names(color_conditions) = levels(metadata$condition)

## ----loaddata, cache=TRUE, results="hide"--------------------------------
this=fcs
# too big to store on Github
db.name = "meta/Spitzer_blood_d3.sqlite3"
this$connectDb(db.name)
this$total.projects=dbListTables(this$conn)

### sort out metadata from database first
# index for metadata
idx=grep("markerIdentity|colnameIndex|fileIdentity|fileIndex|UserHistory|Classification|equipmentInfo|fileComments|SPILL",this$total.projects)

this$df.num = 0
for (i in idx){
  this$df.num = this$df.num + 1 
  this$dataframes.name[this$df.num] = this$total.projects[i]
}

this$fileid_name = "_fileIdentity"
this$markerid_name = "_markerIdentity"

# project vector only
this$total.projects = this$total.projects[-idx]
print(this$total.projects)

### get filenames table
this$current.filetable = this$getDFtable(paste0(this$total.projects[1],this$fileid_name))
this$current.filenames = this$current.filetable[,2]
### get marker names table

this$current.staintable = this$getDFtable(paste0(this$total.projects[1],this$markerid_name))

# file_id same in database and external table
idx.d3 = which(metadata$original_file_name%in%this$current.filenames)
metadata.d3 = metadata[idx.d3,]
metadata.d3$condition =factor(metadata.d3$condition, levels =  
                                c("Bl_CD1_d3","Bl_aPD1_d3","Bl_B6_d3","Bl_untr_d3"))
metadata.d3$sample_id = factor(metadata.d3$sample_id, levels = 
                                 c("aPD1_d3_Bl1","aPD1_d3_Bl2","aPD1_d3_Bl3","B6_d3_Bl1", "B6_d3_Bl2",  "B6_d3_Bl3",  
                                   "CD1_d3-2_Bl1", "CD1_d3_Bl1","CD1_d3_Bl2", "CD1_d3_Bl3", 
                                    "untr_d3_Bl1","untr_d3_Bl2", "untr_d3_Bl3"))

### load matrix data.gated.all
data.gated.all = readRDS(file = sprintf("Rdata/blood_d3_gated_cof%s_bio_trimmed.rds",cofactor))


## ----triplots quadrants, cache=TRUE-----------------------------------------
len.var = ncol(data.gated.all)-1
colvec = colnames(data.gated.all)[2:ncol(data.gated.all)]
len.col = length(colvec)

### initate
registerDoParallel(cl)
it = 0
#quad.df = label.df = data.frame( matrix(
quad.df = data.frame( matrix(
  ,nrow = nrow(metadata.d3)
  ,ncol = 40455
  ,byrow = TRUE
), stringsAsFactors = FALSE
)

ptm <- proc.time()
quad.sample_id = vector()
for ( i in 1:nrow(metadata.d3)) {
  quad.sample_id = c(quad.sample_id,as.character(metadata.d3$sample_id[i]))

  quad.file = vector()
  for ( v1 in 1:(len.col-1) ) {
    it = it +1

    quad.oper <- foreach ( v2=(v1+1):len.col,.combine=cbind )  %dopar% {
      quadrant.vec = vector()

      for ( v3 in 1:len.col ) {
        if ( all(v3 != c(v1,v2)) ) {
          sampl.data = data.gated.all[which(data.gated.all$file_id==metadata.d3$sample_id[i])
                                    ,c(colvec[v1],colvec[v2],colvec[v3])]
          ### NEW::ONLY rows where used if v1 or v2 are >0
          sampl.data = sampl.data[which(sampl.data[,1]>0 & sampl.data[,2]>0),]
          
          # calculate triplot quadrants ---------------------------------------------
          quad.results = this$calc_triplot_quadrant(temp.data = sampl.data, calc.meth = checkCALC, min.cells = min.cells, min.bin.count = bin.min.count)
          quadrant.vec = c(quadrant.vec, quad.results)


          
        }
      }

      return(quadrant.vec)
    }
    
    ### NEW::DONT LOOK AT Q1
    quad.file = c(quad.file,as.vector(quad.oper))
    
    if ( i==nrow(metadata.d3) ) {
      
      label.file = vector()
      for ( v1 in 1:(len.col-1) ) {
        
        label.oper <- foreach ( v2=(v1+1):len.col,.combine=cbind )  %dopar% {
          label.vec = vector()
          
          for ( v3 in 1:len.col ) {
            if ( all(v3 != c(v1,v2)) ) {
              label.vec = c(label.vec,
                            #paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",checkCALC,".","Q1"),
                            paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",checkCALC,".","Q2"),
                            paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",checkCALC,".","Q3"),
                            paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",checkCALC,".","Q4")
              )
            }
          }
          
          label.vec
        }
        label.file = c(label.file,as.vector(label.oper))
      }
    }
    
    printf("%s::%s::quadrants::%s::v1=%s[%s/%s] ready (it=%s)", i, metadata.d3$sample_id[i], checkCALC, colvec[v1],v1,len.col, it)
    print(proc.time() - ptm)
  } 
  quad.df[i,] = quad.file
  
  printf("File %s ready (it=%s)", i, it)
  print(proc.time() - ptm)
}

colnames(quad.df) = label.file
rownames(quad.df) = quad.sample_id

#saveRDS(quad.df, file = sprintf("Rdata/d3_blood_triploT_quadrants_%s_cof%s_bio_cut0.rds",checkCALC,cofactor))
saveRDS(quad.df, file = sprintf("Rdata/d3_blood_triploT_quadrants_%s_cof%s_minbin%s_bio_cut0_trimmed.rds",checkCALC,cofactor,bin.min.count))

print("Done creating blood quadrant table. Look in Rdata/!")

stopCluster(cl)



count = 0
count20 = 0
keep.idx = vector()
for ( i in 1:ncol(quad.df)) {
  keep = TRUE
  if (any(is.na(quad.df[,i]))) {
    printf("col=%s::%s",i,colnames(quad.df)[i])
    count = count + 1
    
    if (sum(is.na(quad.df[,i]))> nrow(quad.df)*0.2) {
      count20 = count20 + 1
      keep = FALSE
    }
  }
  
  if (keep) keep.idx = c(keep.idx,i)
}

printf("keeps=%s out of %s. NAs>20%%=%s",length(keep.idx),ncol(quad.df),count20)
saveRDS(quad.df[,keep.idx], file = sprintf("Rdata/d3_blood_triploT_quadrants_%s_cof%s_minbin%s_keep_bio_cut0_trimmed.rds",
                                           checkCALC,cofactor,bin.min.count))











