#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
rm(list = ls())
source("R/spitzer.r")



this=fcs

cofactor = 0.1

## Metadata
metadata = read.csv(file="meta/Spitzer_blood.csv")
metadata$condition =factor(metadata$condition, levels =  
                             c("Bl_CD1_d3","Bl_aPD1_d3","Bl_B6_d3","Bl_untr_d3"))

color_conditions = c("#006837","#b2182b","#66bd63","#f46d43")
names(color_conditions) = levels(metadata$condition)

### Panel info
panel.tbl = read.csv(file="meta/Panel_for_CyTOF.csv",sep="\t")
protein.names = readRDS("meta/protein_names")

### uninteresting marker
protein.unint = c("CD8","IgM","IgD","CD19","B220","F4-80","FcER1a","PyMT","NK1.1","Ter119")
protein.names = protein.names[-which(protein.names %in% protein.unint)]


## Load CyTOF samples of day3 after treatment
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

### get filenames table
this$current.filetable = this$getDFtable(paste0(this$total.projects[1],this$fileid_name))
this$current.filenames = this$current.filetable[,2]

### get marker names table
this$current.staintable = this$getDFtable(paste0(this$total.projects[1],this$markerid_name))
this$selected.vars = this$getVariables(index=1)
updated.vars = gsub("[^[:alnum:]]","",toupper(this$selected.vars))

protein.tmp = gsub("[^[:alnum:]]","",toupper(protein.names))

if ( all(protein.tmp %in% updated.vars) ) {
  biol.idx = which(updated.vars %in% protein.tmp)
  printf("%s biological marker in total.", length(biol.idx))
} else {
  print("Check metadata staining list!")
}

# file_id same in database and external table
idx.d3 = which(metadata$original_file_name%in%this$current.filenames)
metadata.d3 = metadata[idx.d3,]
metadata.d3$condition =factor(metadata.d3$condition, levels =  
                           c("Bl_CD1_d3","Bl_aPD1_d3","Bl_B6_d3","Bl_untr_d3"))
metadata.d3$sample_id = factor(metadata.d3$sample_id, levels = 
                          c("aPD1_d3_Bl1","aPD1_d3_Bl2","aPD1_d3_Bl3","B6_d3_Bl1", "B6_d3_Bl2",  "B6_d3_Bl3",  
                            "CD1_d3-2_Bl1", "CD1_d3_Bl1","CD1_d3_Bl2", "CD1_d3_Bl3", "untr_d3_Bl1","untr_d3_Bl2", 
                            "untr_d3_Bl3"))

### long version to manipulate data and then merge
if (TRUE) {
  ################################# project: MS Dataset 01
  printf("Load gated data from project: %s", this$total.projects[1])
  
  ### load single files and merge into one data.base.all
  data.gated.all = data.man.all = data.frame()
  ncells.trimmed.total = 0
  for ( i in 1:length(metadata.d3$original_file_name) ) {
  # for ( i in 1:1 ) {
    printf("Processing file %s #%s/%s..",
           metadata.d3$sample_id[i],i,length(metadata.d3$original_file_name))
    db_index = which(this$current.filenames==metadata.d3$original_file_name[i])
    
    ### get marker names and data
    this$selected.vars = this$getVariables(index=db_index)
    data.gated = this$getData(table=this$total.projects[1],fileidx=db_index,cofactor=cofactor)
    
    ### Select biological marker
    updated.vars = gsub("[^[:alnum:]]","",toupper(this$selected.vars))
    if ( all(protein.tmp %in% updated.vars)) {
      biol.idx = which(updated.vars %in% protein.tmp)
      printf("%s biological marker in total.", length(biol.idx))
      data.gated = data.gated[,biol.idx]
    } else {
      print("Check metadata staining list!")
    }
    
    ### need to trim top outliers
    trim.size = 0.0005
    ncells.trimmed = 0
    for (t in 1:ncol(data.gated)) {
      trim.idx = which(data.gated[,t] > quantile(data.gated[,t],c(1-trim.size)))
      data.gated = data.gated[-trim.idx,]
      
      ncells.trimmed = ncells.trimmed + length(trim.idx)
    }
    
    ncells.trimmed.total = ncells.trimmed.total + ncells.trimmed
    printf("Data trimmed: -%s cells", ncells.trimmed)
    
    ### manipulate data part #1
    file_ids = rep(metadata.d3$file_id[i], nrow(data.gated))
    data.man = data.frame(file_id = file_ids, data.gated)
    
    ### concenate data.gated into one matrix: data.gated.all
    # generate sample IDs corresponding to each cell in the 'data' matrix
    col.fileid = rep(metadata.d3$sample_id[i],nrow(data.gated))
    data.gated = cbind(file_id=col.fileid,data.gated)
    if ( i==1 ) {
      data.gated.all = data.gated
    }  else {
      data.gated.all <- merge(data.gated.all,data.gated,all=T)
    }
    
    ### manipulate data part #2
    data.man = melt(data.man, id.var = "file_id", value.name = "expression", variable.name = "marker")
    #head(data.man,2)
    matchmeta = match(data.man$file_id, metadata.d3$file_id)
    data.man$condition = metadata.d3$condition[matchmeta]
    if ( i==1 ) {
      data.man.all = data.man
    }  else {
      data.man.all <- merge(data.man.all,data.man,all=T)
    }
  }
  printf("Total cells trimmed: %s",ncells.trimmed.total)
  print("Done loading and manipulating data for Day3.")
  dim(data.gated.all)
  #[1] 345605     32
  saveRDS(data.gated.all, file = sprintf("Rdata/blood_d3_gated_cof%s_trimmed_bio.rds",cofactor))
  
  dim(data.man.all)
  #[1] 10781428        4
  saveRDS(data.man.all, file = sprintf("Rdata/blood_d3_manipulated_cof%s_trimmed_bio.rds",cofactor))
} else {
  ### load data.gated.all
  data.gated.all = readRDS( file = sprintf("Rdata/blood_d3_gated_cof%s_trimmed_bio.rds",cofactor))
  dim(data.gated.all)
  #[1] 345605     32
  # load matrix data.all
  data.man.all = readRDS(file = sprintf("Rdata/blood_d3_manipulated_cof%s_trimmed_bio.rds",cofactor))
  dim(data.man.all)
  #[1] 10781428        4
}



