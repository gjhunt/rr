# Paper Plot of Single Well
single_well_maker = function(plate,well,add_ecmp_markers=FALSE){
    W = lapply(TRANSFORMATIONS,function(tns)pp_exs$list[[1]]$plates[[tns]]$wells[[plate]][[well]])

    make_inds = function(ecmp){
        ind = data.frame(arrayInd(as.integer(names(fctrs$ecmp$factor)[which(fctrs$ecmp$factor%in%ecmp)]),.dim=c(20,35)))
        #ind = data.frame(arrayInd(c(1,315,320,401,520,700),.dim=c(20,35)))
        colnames(ind) = c("y","x")
        ind$value = 0
        return(ind)
    }

    if(add_ecmp_markers)
        W=lapply(W,function(w)w+geom_point(data=make_inds(c("NID1_1","ELN_3")),mapping=aes(x=y,y=x),color='orange',pch=1,size=7,stroke=2))
    
    W = lapply(seq_along(W),function(i)W[[i]]+
                                       ggtitle(TRANSFORMATIONS[i])+nostrip+
                                       theme(plot.title = element_text(size=35,
                                                                       margin=unit(c(.1,0,.1,0),"cm")),
                                             plot.background = element_rect(fill = "grey")))
    W = lapply(seq_along(W),function(i){
        gb <- rectGrob(height = 1.0, width = 1.0, gp = gpar(lwd = 2, col = "black", fill = NA))
        gt <- gTree(children = gList(ggplotGrob(W[[i]]), gb))
        return(gt)
    })
    names(W) = TRANSFORMATIONS
    return(plot_grid(plotlist=W,nrow=1))
}

p = single_well_maker(5,7)
pdf(paste0(plotdir,"single_well_5_7.pdf"),height=5,width=20)
print(p)
dev.off()

p = single_well_maker("LI8X00515",3)
pdf(paste0(plotdir,"single_well_missing.pdf"),height=5,width=20)
print(p)
dev.off()

p = single_well_maker(9,3)
pdf(paste0(plotdir,"single_well_9_3.pdf"),height=5,width=20)
print(p)
dev.off()

# With Markers

if(sffx==""){
    p = single_well_maker(5,7,TRUE)
    pdf(paste0(plotdir,"single_well_5_7_marker.pdf"),height=5,width=20)
    print(p)
    dev.off()
    p = single_well_maker(9,3,TRUE)
    pdf(paste0(plotdir,"single_well_9_3_marker.pdf"),height=5,width=20)
    print(p)
    dev.off()
}

# Paper Plot of Well Batch

sub_batch = function(j){
    g1 = arrangeGrob(grobs=lapply(c(1:4),function(i)pp_exs$list[[1]]$plates[[j]]$wells[[i]][[3]]+nostrip),nrow=1)
    g2 = arrangeGrob(grobs=lapply(c(9:12),function(i)pp_exs$list[[1]]$plates[[j]]$wells[[i]][[3]]+nostrip),nrow=1)
    g3 = arrangeGrob(grobs=list(textGrob("Batch 1",rot=90,gp = gpar(fontsize = 25)),
                                g1,
                                textGrob("Batch 2",rot=90,gp = gpar(fontsize = 25)),
                                g2
                                ),
                     nrow=2,layout_matrix=matrix(c(1,rep(2,15),3,rep(4,15)),nrow=2,byrow=TRUE))
    gb <- rectGrob(height = 1.0, width = 1.0, gp = gpar(lwd = 4, col = "black", fill = NA))
    gt <- gTree(children = gList(g3, gb))
    return(gt)
}

plts = lapply(c(1,2,4,5),sub_batch)
pdf(paste0(plotdir,"batch_single_well.pdf"),height=8,width=15)
print(plot_grid(plotlist=plts,nrow=2,
          labels=c(padt[c(1,2,4,5)]),
          scale=.97,label_size=25,hjust=-.23,vjust=1.4,align='v'))
dev.off()
