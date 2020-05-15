agg_sv=function(k,q,trans,chir='u',SVDs){
    topk = lapply(SVDs,function(x)x[[trans]][[chir]][,1:k,drop=FALSE])
    all_top_k = Reduce("cbind",topk)
    return(svd_na(all_top_k)$u[,1:q,drop=FALSE])
}

pvd = function(k,q,trans,SVDs){
    U = agg_sv(k,q,trans,"u",SVDs)
    V = agg_sv(k,q,trans,"v",SVDs)
    return(list(u=U,v=V))
}


