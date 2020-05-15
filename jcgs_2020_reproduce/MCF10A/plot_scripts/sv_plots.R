library('data.table')
library('ggplot2')
library('reshape2')
library('cowplot')
library('grid')

source('./plot_fns/sv.R')

flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}

sv_plots = function(var,SVDs){
    cat(var,"\n")

    TNS = c("NT","RR")
    
    cat("Right SVs\n")
    rsv = lapply(TNS,function(trans)right_svs(var,SVDs,trans))
    names(rsv) = TNS

    rsv_ecmp_scatter = sv_meta_scatter(var,SVDs,"right",fctrs$ecmp,kvals=1:4)    
    #rsv_row_scatter = sv_meta_scatter(var,SVDs,"right",fctrs$row,sort_key=TRUE,kvals=1:4)
    #rsv_col_scatter = sv_meta_scatter(var,SVDs,"right",fctrs$col,sort_key=TRUE,kvals=1:4)
    rsv_row_scatter = NULL
    rsv_col_scatter = NULL

    rsv_scatter_by_ecmp = list()
    #for(ecmp_ch in levels(fctrs$ecmp$factor)){
    #    fct = fctrs$ecmp
    #    fct$name = ecmp_ch
    #    levels(fct$factor)[levels(fct$factor)!= ecmp_ch] <- "Other"
    #    fct$factor <- factor(fct$factor,levels=c(ecmp_ch,"Other"),ordered=TRUE)
    #    fct$mtx[fct$mtx != ecmp_ch] <- "Other"
    #    rsv_scatter_by_ecmp[[ecmp_ch]] = sv_meta_scatter(var,SVDs,"right",fct,kvals=1:4)
    #}


    return(list(rsv=rsv,
                rsv_ecmp_scatter=rsv_ecmp_scatter,
                rsv_row_scatter=rsv_row_scatter,
                rsv_col_scatter=rsv_col_scatter,
                rsv_scatter_by_ecmp=rsv_scatter_by_ecmp))
    }

save_sv_plots = function(subdir,SVD_list,uid=NULL){
    sv_plts = lapply(names(SVD_list),sv_plots,SVD_list)
    sv_plts = lapply(sv_plts,flattenlist)
    names(sv_plts) = names(SVD_list)
    for(i in seq_along(sv_plts)){
        any_gg = which(sapply(sv_plts[[i]],function(x)any(class(x)==c("ggplot"))))
        if(length(any_gg)>0)
            sv_plts[[i]][any_gg] <- lapply(sv_plts[[i]][any_gg],ggplotGrob)
        save_plots(sv_plts[[i]],plotdir,subdir,names(sv_plts)[i],uid,height=7)
    }
    plt_nms = names(sv_plts[[1]])
    sapply(plt_nms,function(nm)save_plots(lapply(sv_plts,"[[",nm),plotdir,subdir,nm,uid))
    return(sv_plts)
}

wch_rsv = c("rsv.NT","rsv.RR")

all_plts = save_sv_plots("sv/all_SVD/",SVD_no_center)

ex_plts = save_sv_plots("sv/ex_SVD/",SVD_no_center[ex_ss])

#lapply(names(ex_plts),function(nm)ggsave(paste0(plotdir,"sv/ex_SVD/",nm,"_rsv_well.pdf"),arrangeGrob(grobs=ex_plts[[nm]][wch_rsv],nrow=1),height=15,width=35))
lapply(names(ex_plts),function(nm)ggsave(paste0(plotdir,"sv/ex_SVD/",nm,"_rsv_well.pdf"),plot_grid(plotlist=ex_plts[[nm]][wch_rsv],nrow=1,scale=c(.95,.95)),height=15,width=35))

dapi_plts = save_sv_plots("sv/DAPIASVD/",DAPIASVD_no_center)
lapply(names(dapi_plts),function(nm)ggsave(paste0(plotdir,"sv/DAPIASVD/",nm,"_rsv_well.pdf"),arrangeGrob(grobs=dapi_plts[[nm]][wch_rsv],ncol=1),height=45,width=35))

dapi_morph_plts = save_sv_plots("sv/DAPIASVD_morph/",DAPIASVD_morph)
#lapply(names(dapi_plts),function(nm)ggsave(paste0(plotdir,"sv/DAPIASVD/",nm,"_rsv_well.pdf"),arrangeGrob(grobs=dapi_plts[[nm]][wch_rsv],ncol=1),height=45,width=35))
lapply(names(dapi_plts),function(nm)ggsave(paste0(plotdir,"sv/DAPIASVD/",nm,"_rsv_well.pdf"),plot_grid(plotlist=dapi_plts[[nm]][wch_rsv],nrow=1,scale=c(.95,.95)),height=15,width=35))


if(FALSE){
    ## Mini SVs
    # Mini intensity
    which_ind = which(names(SVD_no_center)=="Cytoplasm_CP_Intensity_IntegratedIntensity_Dapi")
    rsv = lapply(TNS,function(trans)right_svs(which_ind,SVD_no_center,trans))
    rsv = lapply(rsv,function(r){
        gb <- rectGrob(height = 1.0, width = 1.0, gp = gpar(lwd = 4, col = "black", fill = NA))
        gt <- gTree(children = gList(ggplotGrob(r), gb))
    })
    names(rsv) = TNS
    ggsave(paste0(plotdir,"sv/mini_intensity.pdf"),plot_grid(plotlist=rsv),
           width=20,height=8)

    # Mini ASVD
    rsv = lapply(TNS,function(trans)right_svs(1,DAPIASVD_no_center,trans))
    rsv = lapply(rsv,function(r){
        gb <- rectGrob(height = 1.0, width = 1.0, gp = gpar(lwd = 4, col = "black", fill = NA))
        gt <- gTree(children = gList(ggplotGrob(r), gb))
    })
    names(rsv) = TNS
    ggsave(paste0(plotdir,"sv/mini_asvd.pdf"),plot_grid(plotlist=rsv),
           width=20,height=8)
}
