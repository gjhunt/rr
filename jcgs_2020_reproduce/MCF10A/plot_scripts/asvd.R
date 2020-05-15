average_svd_na = function(Y_list,nu=NULL,nv=NULL){
    left_list = lapply(Y_list,gram,project=FALSE)
    right_list = lapply(Y_list,gram,"right",project=FALSE)

    left_avg = mean_mtx_na(left_list)
    left_avg_spsd = spsd_project(left_avg)

    right_avg = mean_mtx_na(right_list)
    right_avg_spsd = spsd_project(right_avg)

    eU = eigen(left_avg_spsd,symmetric=TRUE)
    eV = eigen(right_avg_spsd,symmetric=TRUE)

    if (is.null(c(nu, nv))) 
        nu <- nv <- min(dim(Y_list[[1]]))
    if (is.null(nu)) 
        nu <- sum(eU$values > 0)
    if (is.null(nv)) 
        nv <- sum(eV$values > 0)

    return(list(u = eU$vectors[, 1:nu], du = eU$values[1:nu], v = eV$vectors[, 1:nv], 
                dv = eV$values[1:nv], missing = NULL))
}

stdize_list = function(Y_list){
    sdt_mtxs = function(Yl){
        #sYl = lapply(Yl,function(x)apply(x,2,function(y)(y-mean(y,na.rm=TRUE))))
        sYl = lapply(Yl,function(x)x/sqrt(mean(x^2,na.rm=TRUE)))
        names(sYl) = names(Yl)
        return(sYl)
    }
    Y_list_ret = lapply(Y_list,sdt_mtxs)
    names(Y_list_ret) = names(Y_list)
    return(Y_list_ret)
}
