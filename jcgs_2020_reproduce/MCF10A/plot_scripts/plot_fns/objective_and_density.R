


save_obj_plots = function(scaledY,uid=""){
    
    make_obj_plts = function(s,nm=""){
        cat(nm,"\n")
        T_name = s$T_name
        opt = s$opts[[T_name]]
        l_tried = melt(lapply(opt$opts,"[[","l_tried"))
        obj_tried = melt(lapply(opt$opts,"[[","obj_tried"))
        lh = s$lambda_hat
        obj_df = data.table(Lambda = l_tried$value,
                            OBJ = obj_tried$value,
                            Col = obj_tried$L1)            
        obj_df = obj_df[is.finite(obj_df$OBJ),]
        offset = median(obj_df$OBJ,na.rm=TRUE)
        obj_smooth_df = obj_df[,.(OBJ=offset+OBJ-min(OBJ),Lambda=Lambda),by=.(Col)]

        # Unified Objective
        obj_plt = ggplot(data=cbind(obj_df,z=1),mapping=aes(x=Lambda,y=OBJ,z=z)) + stat_summary_hex(fun=function(z){log(sum(z))}) + scale_fill_gradientn (colours = c("lightgrey","darkblue")) +  
            #geom_smooth(data=cbind(obj_smooth_df,z=1),size=3,inherit.aes=TRUE) +
            geom_vline(xintercept=lh,col="red",size=2)
        obj_plt = obj_plt + labs(title=nm,
                                    subtitle=paste0("log(OBJ) v. Lambda: lambda-hat = ",round(lh,3)," (",T_name,")"),y="log(OBJ)")
        obj_plt = obj_plt+theme(plot.margin = unit(c(1,.5,.5,.5), "cm")) + guides(fill=guide_legend(title="log(count)"))

        # Each Objective Curve
        objs_plts = ggplot(data=obj_df,mapping=aes(x=Lambda,y=OBJ)) + geom_line(data=cbind(obj_df,z=1),mapping=aes(x=Lambda,y=OBJ,group=Col),lty=1,color="black") +  geom_vline(xintercept=lh,col="red",size=2) + 
            labs(title=nm,subtitle=paste0("log(OBJ) v. Lambda: lambda-hat = ",round(lh,3)," (",T_name,")"),y="log(OBJ)")#+ geom_smooth(data=cbind(obj_df,z=1),size=3,inherit.aes=TRUE) 

        # Lambdas Histogram
        ldf = data.table(Lambda = s$lambdas)
        l_plt = ggplot(data=ldf,mapping=aes(x=Lambda))+geom_density(lwd=1)+geom_vline(xintercept=lh,col="red",size=2)
        l_plt = l_plt + labs(title=nm,
                             subtitle=paste0("Density v. Lambda: lambda-hat = ",round(lh,3)," (",T_name,")"))


        # Y v. Lambdas
        lambdas = unname(s$lambdas)
        ds = lapply(TRANSFORMATIONS,function(tns)data.frame(cbind(lambdas,t(s[[tns]]))))
        names(ds) = TRANSFORMATIONS
        mds = lapply(names(ds),function(dn)cbind(melt(ds[[dn]],id.vars="lambdas"),dn))
        D = Reduce("rbind",mds)
        colnames(D) = c("lambdas","variable","value","Transformation")
        Y_l_plt = ggplot(data=cbind(D,z=1),mapping=aes(x=lambdas,y=value,z=z))+ stat_summary_hex(fun=function(z){log(sum(z))}) + scale_fill_gradientn(colours = c("lightgrey","darkblue")) + facet_wrap(~Transformation,scales='free',nrow=1)+geom_vline(xintercept=lh,col="red",size=2) + labs(title=nm,subtitle=paste0("Values v. Lambda: lambda-hat = ",round(lh,3)," (",T_name,")")) + guides(fill=guide_legend(title="log(count)"))

        return(list(obj=obj_plt,objs=objs_plts,l=l_plt,yl=Y_l_plt))
    }

    plts = lapply(seq_along(scaledY),function(i)make_obj_plts(scaledY[[i]],names(scaledY)[i]))

    obj_plts = lapply(plts,"[[","obj")
    save_plots(obj_plts,dir=plotdir,subdir="obj_dens/",name="objective",uid=uid)

    objs_plts = lapply(plts,"[[","objs")
    save_plots(objs_plts,dir=plotdir,subdir="obj_dens/",name="objective_curves",uid=uid)

    lambdas_plts = lapply(plts,"[[","l")
    save_plots(lambdas_plts,dir=plotdir,subdir="obj_dens/",name="lambda_density",uid=uid)

    y_lambdas_plts = lapply(plts,"[[","yl")
    save_plots(y_lambdas_plts,dir=plotdir,subdir="obj_dens/",name="y_lambda_density",uid=uid)
}
