library('scaledown')
library('reshape2')
library('ggplot2')
library('data.table')

Ys = readRDS('../data/Ys')

scaledown_rm = function(N,nme,MELT=FALSE){
    cat(N,"\n")
    mty = as.matrix(Ys[[nme]])
    rm_head = head(order(mty,decreasing=TRUE,na.last=TRUE),n=N)
    rm_tail = head(order(mty,decreasing=FALSE,na.last=TRUE),n=N)
    mty[c(rm_head,rm_tail)] <- NA
    if(MELT)
        mty = melt(mty)[,3,drop=FALSE]
    sdn = scaledown::scaledown(Y=mty,
                               trans_list = list(bcn=scaledown::box_cox_negative),
                               lims_list=list(c(-100,100)),
                               run_parallel=TRUE,verbose=TRUE)
    return(sdn)
}

rm_seq = c(0:5,10,25,50,75,100,500,1000)

sens_compare = function(nme){
    print(nme)
    melted_sdn = lapply(rm_seq,scaledown_rm,nme,TRUE)
    sdn = lapply(rm_seq,scaledown_rm,nme,FALSE)

    melted_lh = sapply(melted_sdn,"[[","lambda_hat")
    lh = sapply(sdn,"[[","lambda_hat")

    d = data.frame(univariate_box_cox=melted_lh,robust_scaledown=lh,N=rm_seq)
    d = melt(d,id.vars='N')
    colnames(d) = c("N","Method","lh")
    levels(d$Method) <- c("Vectorized Lambda","Robust Lambda")
    plt = ggplot(data=d,mapping=aes(x=N,y=lh,color=Method,group=Method))+geom_line(lwd=1)+geom_point()+ggtitle(nme)+labs(x="N",y="Lambda-hat")

    ggsave(paste0('../plots/sens/lambda_',nme,".pdf"),plt)

    TRANS = function (lambda_hat, Y, center = TRUE, scale = TRUE, rm_q = 4, re_estimate = TRUE){
        Yt <- scaledown::box_cox_negative$T(Y, lambda_hat)
        Yt <- as.matrix(Yt)
        if (re_estimate) {
            Yw <- scaledown:::winsor(Yt, 0.001)
            mu <- mean(Yw, na.rm = TRUE)
            Ywc <- Yw - mu
            norm <- sqrt(mean(Ywc^2, na.rm = TRUE))
        }
        else {
            mu <- mu_global
            norm <- norm_global
        }
        if (center) 
            Yt <- Yt - mu
        if (scale) 
            Yt <- Yt/norm
        Yt[abs(Yt) > rm_q] <- NA
        return(Yt)
    }

    make_df = function(tns){        
        if(tns=="Y_bc")
            trans = function(l,x)TRANS(l,x,center=FALSE,scale=FALSE,rm_q=Inf)
        if(tns=="Yt")
            trans = function(l,x)TRANS(l,x)
        if(tns=="Y")
            trans=function(l,x)identity(x)
        if(tns=="Y_bcn")
            trans = function(l,x)TRANS(l,x,rm_q=Inf)
        if(tns=="Y_bcc")
            trans = function(l,x)TRANS(l,x,scale=FALSE,rm_q=Inf)
        mlt=melt(lapply(seq_along(lh),function(i)trans(sdn[[i]]$lambda_hat,sdn[[i]]$Y)))
        ytdf = data.frame(value=mlt$value,run=mlt$L1,
                          Method="Robust")
        mlt=melt(lapply(seq_along(lh),function(i)trans(melted_sdn[[i]]$lambda_hat,melted_sdn[[i]]$Y)))
        melted_ytdf = data.frame(value=mlt$value,run=mlt$L1,
                                 Method="Vectorize")
        ytdf = rbind(ytdf,melted_ytdf)
        ytdf$N = factor(rm_seq[ytdf$run])
        ytdf$lh = lh[ytdf$run]
        ytdf$Method = factor(ytdf$Method)
        ytdf$Trans = tns
        return(ytdf)
    }

    ytdf = rbind(make_df("Y"),make_df("Y_bc"),make_df("Y_bcn"),make_df("Y_bcc"),make_df("Yt"))
    ytdf = data.table(ytdf)
    levels(ytdf$Method) <- c("Robust Lambda","Vectorized Lambda")
    ytdf$Trans = factor(ytdf$Trans)
    levels(ytdf$Trans)=c("No Trans","Box-Cox","Box-Cox Centered","Box-Cox Z-score","Box-Cox Robust Z-score")
    
    dens_plt = ggplot(data=ytdf,mapping=aes(x=value,color=N))+geom_density(aes(y = ..density..),kernel="rectangular",lwd=1)#+geom_histogram(alpha=.25, position="identity",lwd=1)
    dens_plt = dens_plt + facet_wrap(Method~Trans,scales="free",nrow=2)+ggtitle(nme)
    #dens_plt = dens_plt + geom_vline(data=ytdf[,.(max=max(value,na.rm=TRUE)),by=.(Trans,Method,N)],
     #                                mapping=aes(xintercept=max,color=N)) +
     #   geom_vline(data=ytdf[,.(max=min(value,na.rm=TRUE)),by=.(Trans,Method,N)],
     #                                mapping=aes(xintercept=max,color=N))
    
    ggsave(paste0('../plots/sens/dens_',nme,".pdf"),dens_plt,height=20,width=25)

    #qdf = ytdf[,.(q=c(0,.5,1),quant=lapply(c(0,.5,1),function(q)quantile(value,q,na.rm=TRUE))),by=.(N,Trans,Method)]
    #ggplot(data=qdf,mapping=aes(x=as.numeric(N),y=as.numeric(quant),group=interaction(factor(q),Method),color=Method))+geom_line()+facet_wrap(~Trans,scales="free",nrow=1)
    
    
    return(list(d=d,lambda_plt=plt,dens_plt=dens_plt))
}

#run_nme = rev(names(Ys)[which(!sapply(Ys,function(x)any(apply(x,1,function(y)all(is.na(y))))))])
run_nme = names(Ys)
run_nme = run_nme[!sapply(run_nme,function(nme)any(Ys[[nme]]==0,na.rm=TRUE))]
run_nme = run_nme[-(1:15)]
comps = lapply(run_nme,sens_compare)
