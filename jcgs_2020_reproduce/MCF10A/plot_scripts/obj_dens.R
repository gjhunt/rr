# Plate Heatmaps
source('plot_fns/plate_plots.R')
pp_exs = save_plate_plots(scaledY[ex_ss])

lmx = matrix(rep(c(1,rep(2:6,each=10)),2),ncol=2)

dir.create(paste0(plotdir,"obj_dens/plate/"),showWarnings=FALSE)
toss = lapply(names(pp_exs$plot),function(nm)
       ggsave(paste0(plotdir,"obj_dens/plate/",nm,".pdf"),
              arrangeGrob(grobs=c(list(textGrob(nm,gp = gpar(fontsize = 55))),
                                  pp_exs$plot[[nm]]),main=nm,layout_matrix=lmx),height=35,width=34,limitsize=FALSE))


source('single_wells.R')
