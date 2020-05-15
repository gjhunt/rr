library('ggplot2')
source('util.R')
library('grid')
#scaledY=readRDS('../scaledY/scaledY.rds')
#plotdir='../plots/'

make_simple_histogram = function(nme){
    df = data.frame(x=as.vector(scaledY[[nme]]$NT))
    plt = ggplot(data=df,mapping=aes(x=x))+geom_histogram() + labs(x="Value",y="Count")
    plt = plt + ggtitle(paste("Histogram of",nme))
    return(plt)
}


plts = lapply(names(scaledY),make_simple_histogram)
names(plts) = names(scaledY)

save_plots(plts,dir=plotdir,subdir="simple_hist/",name="shist",
           height=5,width=10)
