trans = TRANSFORMATIONS

get_SVD = function(sY,center_fn=identity){
    ret = lapply(trans,function(x)sY[[x]])
    names(ret) = trans
    ret_adj = lapply(ret,center_fn)
    ret_adj = lapply(ret_adj,svd_na)
    return(ret_adj)
}

get_Ys = function(sY,center_fn=identity){
    ret = lapply(trans,function(x)sY[[x]])
    names(ret) = trans
    ret_adj = lapply(ret,center_fn)
    return(ret_adj)
}

get_ASVD=function(Y_list){
    sep = lapply(trans,function(x)
        lapply(Y_list,"[[",x))
    names(sep) = trans
    ret_adj = lapply(sep,average_svd_na)
    return(ret_adj)               
}

get_PVD=function(SVDs,p,q){
    trans = names(SVDs[[1]])
    ret_adj=lapply(trans,function(x)
        pvd(p,q,x,SVDs))
    names(ret_adj) = trans
    return(ret_adj)
}

calc_ASVD=function(Y_list,ss=NULL){
    if(is.null(ss)) ss = 1:length(Y_list)
    ASVD=get_ASVD(Y_list[ss])
    ASVD=list(ASVD=ASVD)
    return(ASVD)
}

calc_PVD=function(SVD_list,ss=NULL,p,q){
    if(is.null(ss)) ss = 1:length(SVD_list)
    PVD=get_PVD(SVD_list[ss],p,q)
    PVD=list(PVD=PVD)
    return(PVD)
}
