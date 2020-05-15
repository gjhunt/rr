library('ggplot2')

theme_set(theme_bw())

library('R.utils')
#sffxes = c("","_no_nid_eln","_no_nid_eln_thbs_fn")
args = commandArgs(asValues=TRUE)
if(is.null(args$sffx)){
    sffx <- ""
} else {
    sffx <- args$sffx
}
cat("sffx:",sffx,"\n")

source('util.R')
source('./read_sy.R')
ex_ss = c("Cells_CP_AreaShape_Area","Cells_CP_AreaShape_Compactness","Spot_PA_SpotCellCount","Cytoplasm_CP_Intensity_IntegratedIntensity_Dapi")

source('simple_scatter.R')

### Objective and Density Plots

source('density_plot.R')
source('obj_dens.R')

### Canonical Correlation Plots
source('svd_na.R')
source('pvd.R')
source('asvd.R')
source('get_util.R')

# Ys
Y_no_center=lapply(scaledY,get_Ys,identity)

# Individual
SVD_no_center=lapply(scaledY,get_SVD,identity)

# DAPI
dapi_ss = which(sapply(lapply(Y_no_center,"[[","NT"),function(x)!any(apply(x,1,function(y)all(is.na(y))))))
DAPIASVD_no_center=calc_ASVD(Y_no_center,ss=dapi_ss)

dapi_levs = lapply(lapply(Y_no_center[dapi_ss],"[[","NT"),svd_na)
dapi_levs = sapply(seq_along(dapi_levs),function(i)max(diag(dapi_levs[[i]]$v[,1:10]%*%t(dapi_levs[[i]]$v[,1:10]))))
dapi_hl_ss = dapi_ss[order(dapi_levs,decreasing=TRUE)[1:5]]
DAPIASVD_hl=calc_ASVD(Y_no_center,ss=dapi_hl_ss)

dapi_morph_ss = dapi_ss[grep("AreaShape",names(dapi_ss))]
DAPIASVD_morph=calc_ASVD(Y_no_center,ss=dapi_morph_ss)

#DAPIPVD_no_center=calc_PVD(SVD_no_center,ss=dapi_ss,p=192,q=192)
source('./cca_plots.R')

## Singular Vector Plots
source('./sv_plots.R')
