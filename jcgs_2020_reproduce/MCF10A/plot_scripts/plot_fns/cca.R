library('ggplot2')
library('grid')
library('data.table')
library('cowplot')
library('CCA')
library('MASS')
library('reshape2')

make_cca = function(batch,SVDs,add_SVDs=NULL,chir="u",ncor=NULL,kmax=NULL){
    cat(batch$name,"\n")
    batch_basis = batch$basis

    denom = sum(svd(scale(batch_basis,scale=FALSE))$d>1E-10)
    if(is.null(kmax))
        kmax = ncol(SVDs[[1]][[1]][[chir]])

    ccaout = function(k,SVDs,trans,var){
        W = SVDs[[var]][[trans]][[chir]][,1:k,drop=FALSE]
        mss = SVDs[[var]][[trans]]$missing
        batch_mtx = batch$mtx
        batch_mtx[mss] <- NA
        lvs = unique(as.vector(batch_mtx))
        lvs = lvs[!is.na(lvs)]
        keep_rows = which(!apply(batch_mtx,1,function(x)all(is.na(x))))
        if(length(lvs)==1)
            return(NA)
        W = W[keep_rows,,drop=FALSE]
        batch_basis = batch_basis[keep_rows,]
        denom = sum(svd(scale(batch_basis,scale=FALSE))$d>1E-10)
        C = tryCatch({
            cancor(W,batch_basis)
        },error=function(e){
            list(cor=NA)
        })
        denom=min(denom,ncor)
        if(is.null(ncor))
            ncor=length(C$cor)
        ncor = min(ncor, k)
        return(sum(C$cor[1:ncor]^2)/denom)
    }

    synames=names(SVDs)
    SVD_out = function(SVDs){
        ot = lapply(seq_along(SVDs),function(i){
            cat(i,"\n")
            km = min(kmax,ncol(SVDs[[1]][[1]][[chir]]))
            L = lapply(names(SVDs[[i]]),function(trans)sapply(1:(km),ccaout,SVDs,trans,i))
            names(L) = names(SVDs[[i]])
            M = cbind(melt(data.frame(L,check.names=FALSE),id.vars=NULL),1:(km))
            colnames(M) = c("Transformation","value","index")
            return(M)
            })
        return(ot)
    }

    cca_out = SVD_out(SVDs)
    names(cca_out) = names(SVDs)

    add_out=NULL
    if(!is.null(add_SVDs)){
        add_out = SVD_out(add_SVDs)
    }

    return(list(cca=cca_out,add=add_out))
}

save_cca_plots = function(cca_out,add_out,batch,SVDs,add_SVDs=NULL,chir="u",ncor=NULL,kmax=NULL,uid=NULL,subdir=NULL,save=TRUE){

    cat(batch$name,"\n")
    batch_basis = batch$basis
        synames=names(SVDs)

    denom = sum(svd(scale(batch_basis,scale=FALSE))$d>1E-10)
    if(is.null(kmax))
        kmax = ncol(SVDs[[1]][[1]][[chir]])


    make_cca_plot = function(i){
        max_ks = sapply(names(SVDs[[i]]),function(x)sum(SVDs[[i]][[x]]$du/SVDs[[i]][[x]]$du[1] > 1E-10))
        max_ks = melt(max_ks)
        max_ks$Transformation = rownames(max_ks)
        mss = SVDs[[i]][[1]]$missing
        batch_mtx = batch$mtx
        batch_mtx[mss] <- NA
        kmax_mss = sum(apply(batch_mtx,1,function(x)!all(is.na(x))))
        max_ks = max_ks[max_ks$value != kmax_mss,]

        cca_out[[i]]$value <- as.numeric(cca_out[[i]]$value)
        cca_out[[i]]$Transformation <- factor(cca_out[[i]]$Transformation,levels=TRANSFORMATIONS,ordered=TRUE)

        if(length(unique(cca_out[[i]]$value))==1)
            return(NULL)


        DFi = data.table(cca_out[[i]])
        
        thinned <- floor(seq(from=1,to=max(DFi$index),length=10))
        
        p = ggplot(data=DFi,mapping=aes(x=index,y=value,color=Transformation,order=Transformation,shape=Transformation))
        p = p + geom_abline(slope=1/kmax_mss,intercept=0,lty=2) + geom_hline(yintercept=1)+geom_vline(xintercept=kmax_mss,lwd=1)
        p = p + geom_vline(data=max_ks,mapping=aes(xintercept=value,color=Transformation))
        p = p + geom_point(data=DFi[index%in%thinned,],size=2.5,stroke=1.5)
        p = p + geom_line(lwd=1,position=position_jitter(w=0, h=0.001))
        p = p + ggtitle(paste0(synames[i],": Mean of Squared Can. Cors. between Batch and First k PCs"))+labs(x="Num. PCs (k)",y="Mean of Squared Canonical Cors.")+ylim(0,1.01)+scale_x_continuous("",limits=c(1,kmax),breaks=scales::pretty_breaks(n = 10))
        if(!is.null(add_out)){
          p = p + geom_line(data=add_out[[i]],lty=2,lwd=2)
        }
        color_map = ggplotColors(length(TRANSFORMATIONS))
        names(color_map) = TRANSFORMATIONS
        p = p + scale_color_manual(breaks=TRANSFORMATIONS,values=color_map)+scale_shape_manual(values=0:10)
        return(p)
    }

    make_auc_plot = function(top=kmax){

        calc_auc=function(i,top=kmax,output){
            dt=data.table(output[[i]])
            dt=dt[,.(auc=mean(head(value,top))),by=.(Transformation)]
            dt$var = synames[i]
            colnames(dt)=c("Transformation","AUC","Var")
            return(dt)
        }

        aucs = lapply(1:length(cca_out),calc_auc,top=top,cca_out)
        auc_df = Reduce(rbind,aucs)
        auc_df = auc_df[complete.cases(auc_df),]
        auc_df$Var = factor(auc_df$Var,ordered=TRUE,
                            levels=unique(auc_df$Var))

        ordering = auc_df[,.(diff=AUC[Transformation==TRANSFORMATIONS[5]]-AUC[Transformation==TRANSFORMATIONS[1]]),by=.(Var)]
        var_ordered = ordering$Var[order(ordering$diff,decreasing=TRUE)]
        auc_df$Var = ordered(auc_df$Var,levels=var_ordered)
        

        if(!is.null(add_out)){
            add_aucs = lapply(1:length(add_out),calc_auc,top=top,add_out)
            add_auc_df = Reduce(rbind,add_aucs)
            add_auc_df[Var %in% auc_df$Var,]
            add_auc_df$Var = ordered(add_auc_df$Var,
                                     levels=levels(auc_df$Var))
            add_auc_df = add_auc_df[complete.cases(add_auc_df),]
            lvs = levels(add_auc_df$Var)
            lvs = gsub("_CP_AreaShape_","...",lvs)
            lvs = gsub("_CP_Intensity_","...",lvs)
            lvs = gsub("Spot_PA_","...",lvs)
            lvs = gsub("Cells_","Cell_",lvs)
            levels(add_auc_df$Var) <- lvs
        }

        auc_df$Transformation <- ordered(auc_df$Transformation,levels=TRANSFORMATIONS)

        lvs = levels(auc_df$Var)
        lvs = gsub("_CP_AreaShape_","...",lvs)
        lvs = gsub("_CP_Intensity_","...",lvs)
        lvs = gsub("Spot_PA_","...",lvs)
        lvs = gsub("Cells_","Cell_",lvs)
        levels(auc_df$Var) <- lvs

        p=ggplot(data=auc_df,mapping=aes(y=AUC,x=Var,color=Transformation,shape=Transformation))+geom_point(size=3,stroke=1)
        if(!is.null(add_out)){
            p=ggplot(data=auc_df,mapping=aes(x=Var,y=AUC,color=Transformation,shape=Transformation))+
                geom_point(data=add_auc_df,mapping=aes(x=Var,y=AUC,shape=Transformation),color='black',size=5)
            +
                geom_line(data=add_auc_df,mapping=aes(x=Var,y=AUC,group=Transformation),color='black',lty=2)+
                geom_point(size=5,stroke=2)
        }
        p=p+theme(axis.text.x = element_text(angle = 75, hjust = 1,size=10))
        p=p+ggtitle(paste0("AUC by Feature (over first ",top,")"))+scale_shape_manual(values=0:10)
        p =p + expand_limits(y=1)+geom_hline(yintercept=.5)+geom_hline(yintercept=1)
        p = p + labs(x="Feature")
        p = p + scale_color_manual(breaks=TRANSFORMATIONS,values=ggplotColors(length(TRANSFORMATIONS)))
        return(p)
    }

    make_auc_plots = function(){
        return(list(all=make_auc_plot(),
                    top=make_auc_plot(floor(kmax/3)),
                    batch=make_auc_plot(denom)))
    }

    cca_plts=NULL
    auc_plts=NULL
    if(save){
        cca_plts = lapply(1:length(cca_out),make_cca_plot)
        names(cca_plts) = names(SVDs)
        cca_plts = cca_plts[!sapply(cca_plts,is.null)]
        auc_plts = make_auc_plots()
        save_plots(auc_plts,plotdir,subdir=paste0("cca/",subdir),"auc_",uid=uid,height=7,width=12)
        save_plots(cca_plts,plotdir,subdir=paste0("cca/",subdir),"cca_",uid=uid)
    }
    return(list(data=cca_out,add_data=add_out,cca_plts=cca_plts,auc_plts=auc_plts))
}


meta_cca_plot=function(exs,SVDs,btch="stain"){

    max_k_feature = function(i){
        L=lapply(TRANSFORMATIONS,function(x)sum(SVDs[[i]][[x]]$du/SVDs[[i]][[x]]$du[1] > 1E-10))
        names(L) = TRANSFORMATIONS
        return(L)
    }
    max_ks = lapply(names(exs),max_k_feature)
    names(max_ks) = names(exs)
    max_ks = melt(max_ks)
    colnames(max_ks) =c("index","Trans","Feature")
    max_ks = data.table(max_ks)
    hard_max = max_ks[,.(index=max(index)),by=.(Feature)]

    d = lapply(seq_along(exs),function(i)data.table(cbind(exs[[i]],names(exs)[i])))
    d = lapply(seq_along(d),function(i)d[[i]][index<hard_max$index[i],])
    D = Reduce("rbind",d)
    colnames(D) = c("Trans","value","index","Feature")

    D$Feature = gsub("CP_","",D$Feature)
    D$Feature = gsub("AreaShape_","",D$Feature)
    D$Feature = gsub("_Intensity_","_",D$Feature)
    D$Feature = gsub("Cells_","Cell_",D$Feature)
    
    max_ks$Feature = gsub("CP_","",max_ks$Feature)
    max_ks$Feature = gsub("AreaShape_","",max_ks$Feature)
    max_ks$Feature = gsub("_Intensity_","_",max_ks$Feature)
    max_ks$Feature = gsub("Cells_","Cell_",max_ks$Feature)

    hard_max$Feature = gsub("CP_","",hard_max$Feature)
    hard_max$Feature = gsub("AreaShape_","",hard_max$Feature)
    hard_max$Feature = gsub("_Intensity_","_",hard_max$Feature)
    hard_max$Feature = gsub("Cells_","Cell_",hard_max$Feature)
    

    thinned <- floor(seq(from=1,to=max(D$index),length=10))
    p=ggplot(data=D,mapping=aes(x=index,y=value,color=Trans,shape=Trans))+geom_point(data=D[index%in%thinned,],size=2)+geom_line(lwd=1,position=position_jitter(w=0, h=0.001))
    p = p+facet_wrap(~Feature,scales="free_x",nrow=1)+geom_vline(data=max_ks,mapping=aes(xintercept=index,color=Trans))+geom_vline(data=hard_max,mapping=aes(xintercept=index))+geom_abline(data=hard_max,mapping=aes(slope=1/index,intercept=0),lty=2) + geom_hline(yintercept=1)+labs(x="Num. Components (k)",y="Mean of Squared Canonical Cors.")+theme(legend.position="bottom")+ scale_color_manual(breaks=TRANSFORMATIONS,values=ggplotColors(length(TRANSFORMATIONS)))
    p = p +scale_shape_manual(values=0:10)
    p

    return(p)
}
