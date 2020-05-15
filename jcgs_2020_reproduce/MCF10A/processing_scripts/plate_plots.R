library('ggplot2')
library('grid')
library('dplyr')
library('data.table')

make_plate_plot = function(var){
    if(is.factor(d[[var]]) & length(levels(d[[var]]))==0)
       return(NULL)
    plt = ggplot(data=d,mapping=aes_string(x="PrintColumn",y="PrintRow",fill=var))+geom_tile()+scale_y_reverse()+facet_grid(WellIndex~Barcode) + ggtitle(var)
    if(is.factor(d[[var]]) & length(unique(d[[var]]))>20)
        plt = plt + theme(legend.position="top")
    return(plt+theme(plot.margin = unit(c(1,1,1,1), "cm")))
}

ss=c(30)
plate_plots = lapply(colnames(d)[ss],make_plate_plot)
names(plate_plots)=colnames(d)[ss]

dir.create('../plate_design/',showWarnings=FALSE)
pdf("../plate_design/ligand_design.pdf",width=20,height=12)
for(p in names(plate_plots)){
    cat(p,"\n")
    grid.draw(plate_plots[[p]])
}
dev.off()
