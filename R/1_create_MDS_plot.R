#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
rm(list = ls())
source("R/spitzer.r")



### need to be created from 0_create_dataframe.R
data.man.all = readRDS(file = sprintf("Rdata/blood_d3_manipulated_cof%s_trimmed_bio.rds",cofactor))

# Get the median marker expression per sample
expr_median_sample_tbl = data.man.all[,c(1,2,3)] %>% group_by(file_id,marker) %>% summarize_all(funs(median))
expr_median_sample_tbl = as.data.frame(expr_median_sample_tbl)
# and create matrix with markers in rows and samples in columns
median.mat = matrix(
  expr_median_sample_tbl[,3]
  , nrow = length(unique(expr_median_sample_tbl$marker))
  , ncol = length(metadata.d3$sample_id)
)
colnames(median.mat) = metadata.d3$sample_id
rownames(median.mat) = unique(expr_median_sample_tbl$marker)


# mds plot samples --------------------------------------------------------
mds = plotMDS(median.mat,plot=F)

ggdf.mds <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                       sample_id = colnames(median.mat))
matchmeta = match(ggdf.mds$sample_id, metadata.d3$sample_id)
ggdf.mds$condition <- metadata.d3$condition[matchmeta]

ggmds = ggplot(ggdf.mds, aes(x = MDS1, y = MDS2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  theme(text = element_text(size=14)) +
  scale_color_manual(values = color_conditions) +
  coord_fixed()
ggsave(sprintf("figure/blood_d3_manipulated_cof%s_bio_mds_samples.pdf",cofactor),
       width=9,height=5)

