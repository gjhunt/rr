# Density Plots
source('plot_fns/density.R')
for(fct in c("well","ligand","stain","plate")){
    toss = save_dens_plots(fctr=fct,scaledY=scaledY[ex_ss],uid=paste0("_",fct))
}
for(fct in c("well","ligand","stain","plate")){
    toss = save_dens_plots(fctr=fct,scaledY=scaledY[ex_ss],uid=paste0("_",fct),
                           suppl=TRUE)
}
  
