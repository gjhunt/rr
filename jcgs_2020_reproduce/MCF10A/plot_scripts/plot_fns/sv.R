library('RColorBrewer')

getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

left_sv_k = function(k,var,SVDs,trans,fact,leg=FALSE){
    fname = fact$name
    f = fact$factor
    svdf = data.frame(sv=SVDs[[var]][[trans]]$u[,k])
    svdf$m = 1:nrow(svdf)
    svdf$f = f
    M = fact$mtx
    M[SVDs[[var]][[TRANSFORMATIONS[1]]]$missing] <- NA
    keep = apply(M,1,function(x)!all(is.na(x)))
    svdf = svdf[keep,]
    plt=ggplot(data=svdf,mapping=aes(x=m,y=sv,color=f,shape=f))+geom_point()+scale_shape_manual(values=c(0:25,65:90,97:122),name=fname,labels=levels(f))+scale_color_manual(name=fname,labels=levels(f),values=getPalette(length(unique(f))))
    if(!leg)
        plt = plt + theme(legend.position="none")
    return(plt)
}

left_svs = function(var,SVDs,trans,fact){
    fname = fact$name
    f = fact$factor
    lgnd = g_legend(left_sv_k(1,var,SVDs,trans,fact,leg=TRUE))
    plt = do.call(plot_grid,c(lapply(1:4,left_sv_k,var,SVDs,trans,fact),list(ncol=2,labels=paste0("Left Singular Vector ",1:4))))
    plt = plot_grid(plt,lgnd,rel_widths=c(1.,.3))
    plt=meta_title(plt,tools::toTitleCase(paste0(var," Left singular values, color by ",fname,": ",trans)),size=25)
    return(plt)
}

right_sv_k = function(k,var,SVDs,trans){
    svdf = data.frame(sv=SVDs[[var]][[trans]]$v[,k])
    svdf$row = apply(fctrs$row$mtx,2,unique)
    svdf$col = apply(fctrs$col$mtx,2,unique)
    plt = ggplot(data=svdf,mapping=aes(x=col,y=row,fill=sv))+geom_tile()+theme(aspect.ratio=1)+scale_fill_gradient2()+scale_y_reverse(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+guides(fill=FALSE)
    plt = plt + nostrip
    return(plt)
}


right_svs = function(var,SVDs,trans){
    wch=1:10
    plt=do.call(plot_grid,c(lapply(wch,right_sv_k,var,SVDs,trans),list(nrow=2)))
    plt=meta_title(plt,paste0(trans,": ",var),size=35)
    return(plt)
}

meta_title = function(p,title,size=14){
    title <- ggdraw() + draw_label(title, fontface='bold',size=size)
    plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
}


sv_meta_scatter = function(var,SVDs,chir,fact,sort_key=FALSE,kvals=1:2,TRANSFORMATIONS=c("NT","RR")){
    
    L = lapply(1:(length(kvals)/2),function(i)sv_scatter(var,SVDs,chir,fact,sort_key=FALSE,kvals=c(2*i-1,2*i),TRANSFORMATIONS=c("NT","RR")))

    aL = lapply(L,function(l)l+theme(legend.position="none"))
    gout = do.call(grid.arrange,aL)
    aL[[length(aL)+1]] <- g_legend(L[[1]])

    plt = grid.arrange(grobs=aL,layout_matrix=t(array(c(1,1,2,2,1,1,2,2,1,1,2,2,3,3,3,3),c(4,4))))

    return(plt)
}


sv_scatter = function(var,SVDs,chir,fact,sort_key=FALSE,kvals=1:2,TRANSFORMATIONS=c("NT","RR")){
    fname = fact$name
    f = fact$factor
    chiruv=switch(chir,"right"="v","u")

    M = fact$mtx
    M[SVDs[[var]][[TRANSFORMATIONS[1]]]$missing] <- NA
    keep = apply(M,switch(chiruv,"v"=2,1),function(x)!all(is.na(x)))

    Ws = lapply(TRANSFORMATIONS,function(trans){
        SV = SVDs[[var]][[trans]][[chiruv]]
        WY = data.frame(SV[,kvals])
        WY$f = f
        WY$trans = trans
        WY = WY[keep,]
        return(WY)
    })
   
    W = Reduce("rbind",Ws)
    W$trans = ordered(W$trans,levels=TRANSFORMATIONS)

    shapevals=c(65,66,14,121,76,6,4,19,108,71,68,23,9,97,5,72,110,18,98,17,2,118,16,25,114,116,20,81,67,1,111,87,13,106,74,100,84,120,88,89,73,86,117,15,105,109,99,12,10,119,85,103,80,70,3,112,75,78,22,107,7,0,82,79,90,104,115,101,8,113,83,102,11,24,21,69,77,122)#c(0:25,65:90,97:122)

    #cat("c(",paste("\"",getPalette(46)[sample(46)],"\"",collapse=",",sep=""),")")
    colorvals = c("#BC6334","#78A619","#1B9E77","#916B80","#9F6867","#786A51","#AFA80D","#D39B09","#937131","#DCA305","#9D704C","#A551A1","#66A61E","#E7298A","#D6338F","#C63D95","#AE664D","#8A6F3C","#C24C6B","#6F685B","#8AA716","#955BA7","#AF7D19","#877A34","#E6AB02","#D43A7A","#6C8344","#B6479B","#D95F02","#AF5E5B","#369566","#CA940D","#A27123","#9CA812","#816C46","#7570B3","#78942D","#A6761D","#CA611B","#8565AD","#C1A909","#666666","#B88515","#D3AA05","#836D99","#C18C11","#9C7327","#518C55","#BD6812","#8A823C")

    if(sort_key){
        shapevals = sort(shapevals)
        colorvals = getPalette(length(colorvals))
    }
    
    
    plt=ggplot(data=W,mapping=aes(x=X1,y=X2,color=f,shape=f))+geom_point(size=5)+geom_point(data=W[W$f!="Other",],size=5)
    svar = var
    svar = rev(strsplit(svar,"_")[[1]])[1]
    plt = plt +facet_wrap(~trans,scales='free')+ggtitle(paste0("RSV ",kvals[1],", ",kvals[2],": ",svar))+scale_shape_manual(values=shapevals,name=fname,labels=levels(f),guide = guide_legend(nrow=6))+scale_color_manual(name=fname,labels=levels(f),values=colorvals)
    plt = plt + theme(strip.text=element_text(size=20),
                      legend.text=element_text(size=15),
                      legend.title=element_text(size=20),
                      legend.position="bottom",
                      title=element_text(size=20)) + 
        labs(x=paste0("SV",min(kvals)),y=paste0("SV",max(kvals)))
    return(plt)
}
