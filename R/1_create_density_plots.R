#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
rm(list = ls())
source("R/spitzer.r")



### need to be created from 0_create_dataframe.R
data.man.all = readRDS(file = sprintf("Rdata/blood_d3_manipulated_cof%s_trimmed_bio.rds",cofactor))

# plot all densitites -----------------------------------------------------
ggp = ggplot(data.man.all, aes(x = expression, color = condition,
                               group = file_id)) +
  geom_density() +
  facet_wrap( ~ marker, nrow = 10, scales = "free") +
  theme_bw() +
  theme( text = element_text(size = 16)) +      
  scale_color_manual(values = color_conditions)

### get plot dimension
n_panels <- length(unique(ggplot_build(ggp)$data[[1]]$PANEL))
plot.dims = wrap_dims(n_panels)
### save boxplot as pdf
ggsave(sprintf("figure/blood_d3_manipulated_cof%s_bio_dens.pdf",cofactor),
       width=3*plot.dims[1],height=3*plot.dims[2])
