library('ggplot2')
library('reshape2')
library('data.table')
library('ggpubr')
library('grid')
library('gridExtra')
library('cowplot')
library('viridis')

TRANSFORMATIONS = c("NT","G","Z","O","RR")
padt = c("NT"," G"," Z"," O","RR")

## Set Plotting Directory
plotdir = paste0("../plots",sffx,"/")
dir.create(plotdir,showWarnings=FALSE)

## Read in Scaled Data and throw out optimization results
syd = paste0("../scaledY/scaledY",sffx,".rds")
scaledY = readRDS(syd)
uniques = sapply(scaledY,function(x)length(unique(unlist(x$NT))))
scaledY = scaledY[uniques > 10]
scaledY = scaledY[!grepl("Gated",names(scaledY))]

## Cannot be estimated
aY = lapply(scaledY,"[",TRANSFORMATIONS)
unestimatable = which(sapply(aY,function(x)any(sapply(x,function(y)all(!is.finite(y))))))
rm(aY);gc();
cat("Cannot use these features:\n")
cat(paste(names(unestimatable),"\n"))
scaledY = scaledY[-unestimatable]
rm_opts = function(x){
    return(x[!(names(x)%in%"opts")])
}
scaledY = lapply(scaledY,rm_opts)
gc();

## Read in Factors
fctrs = readRDS(paste0("../data/fctrs",sffx,".rds"))
