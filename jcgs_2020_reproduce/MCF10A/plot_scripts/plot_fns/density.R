save_dens_plots = function(fctr,scaledY,uid="",suppl=FALSE){

    if(suppl)
        uid=paste0(uid,"_suppl")

    batch_factor = fctrs[[fctr]]$factor

    make_dens_plts = function(sdn,nm=""){
        cat(nm,"\n")
        ## Density Plots
        y_batch = function(Ytmp,bf=batch_factor){
            Ytmp = cbind(bf,Ytmp)
            Ytmp = data.frame(Ytmp)
            Ytmp = melt(Ytmp,id='bf')
            y_df  = data.table(batch=Ytmp$bf,Y=Ytmp$value)
            return(y_df)
        }

        trans = TRANSFORMATIONS
        y_batches = lapply(trans,function(x)y_batch(sdn[[x]]))
        y_batches = do.call("cbind",lapply(y_batches,"[[","Y"))
        colnames(y_batches) = trans
        y_batches = data.table(y_batches)
        y_batches$batch = rep(batch_factor,nrow(y_batches)/length(batch_factor))
        
        y_df = melt(y_batches,id="batch",variable.name="Trans")
        y_df$batch = factor(y_df$batch)
        levels(y_df$batch) = levels(batch_factor)
        lh = sdn$lambda_hat
        labs <- TRANSFORMATIONS
        #labs[c(2,5)] = paste(labs[c(2,5)],"\n",sdn$T_name,"(",round(lh,3),")",sep="")
        levels(y_df$Trans) <- labs

        NROW = 1
        if(suppl)
            NROW=1
        
        d_plt = ggplot(data=y_df,mapping=aes(x=value))+geom_density(lwd=2) + facet_wrap(~Trans,scales='free',nrow=NROW)
        d_plt = d_plt + geom_density(data=y_df,mapping=aes(x=value,group=batch,fill=batch,lty=batch),alpha=.25)+guides(fill=guide_legend(title=tools::toTitleCase(fctrs[[fctr]]$name)),lty=FALSE)
        d_plt = d_plt + theme(legend.position="bottom",
                              axis.text.x = element_text(angle=60, hjust=1),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank(),
                              plot.margin=unit(c(0,.1,.1,0),"cm"),
                              panel.spacing = unit(0, "lines"),
                              strip.text.x=element_text(size=20)
                              )+
            labs(title=nm,x=NULL,y=NULL)+
            scale_fill_discrete(breaks=unique(fctrs[[fctr]]$factor),drop=!suppl)
        leg=NULL
        if(suppl){
            leg = g_legend(d_plt)
            d_plt = d_plt + theme(legend.position="right",
                                  axis.text.x = element_text(angle=90, hjust=1),)+
                guides(fill=FALSE)
        }
        return(list(dens=d_plt,leg=leg))
    }
    
    plts = lapply(seq_along(scaledY),function(i)make_dens_plts(scaledY[[i]],names(scaledY)[i]))
    dens_plts = lapply(plts,"[[","dens")
    if(suppl){
        dens_plts = lapply(dens_plts,ggplotGrob)
        p = arrangeGrob(grobs=c(dens_plts,list(plts[[1]]$leg)),ncol=1)
        pdf(paste0(plotdir,"obj_dens/dens",uid,".pdf"),height=12,width=10)
        grid.draw(p)
        dev.off()
    } else {
        save_plots(dens_plts,dir=plotdir,subdir="obj_dens/",name="dens",uid=uid,
                   height=4,width=10)
    }
    return(dens_plts)
}

