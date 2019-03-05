#' Re-scale a data matrix
#' 
#' This transformation is three steps (1) Gaussianize the data, (2) z-score Transform the data, and (3) remove extreme outliers from the data. The sequence of these transformations helps focus further analyses on consequential variance in the data rather than having it be focused on variation resulting from the feature's measurement scale or outliers. 
#' @export
#' @param Y Data matrix to be transformed.
#' @param trans_list List of transformations to be considered. See function list_transformations. Each element of the list should be a list containing the transformation function as the first element and the derivative of the transformation function as the second argument. The first argument of each function should be the data, the second the transformation parameter.
#' @param lims_list List of optimization limits for each transformation from trans_list. This should be a list with one element per transforamtion parmeter. The element of the list for each transforamtion family should be a list of two-element vectors control the limits for each parameter of the transformation family. 
#' @param opt_control Optional optimization controlling parameters for DEoptim control argument. See the DEoptim package for details. 
#' @param z The O-step cutoff value. Points are removed if their robust z-score is above z in magnitude. 
#' @param q The Z-step winsorizing quantile cutoff. The quantile at which to winsorize the data when calculating the robust z-scores. 
#' @param ncores Number of cores to use if running parallel. 
#' @param run_parallel a boolean, if TRUE the method will be run run column-wise in parallel. 
#' @param verbose a boolean, if TRUE then save optimization output in local directory '.rrscale'
#' @param zeros How to deal with zeros in the data set. If set to FALSE the algorithm will fail if it encounters a zero. If set to a number or 'NA' then the zeros are replaced by this number or 'NA'.
#' @param opts Boolean determining if optimization output is returned. Defaults to FALSE. 
#' @param seed Sets the seed before running any other analyses. 
#' @return A list of output:
#' \itemize{
#' \item{opts:} the optimization output for all transformation families and all columns
#' \item{pars:} the optimal parameters for each column for the optimal family
#' \item{par_hat:} the estimated optimal paramter
#' \item{NT:} the original data
#' \item{RR:} the robust-rescaled data
#' \item{G:} gaussianized data
#' \item{Z:} robust z-transformed data
#' \item{O:} data with outliers removed
#' \item{rr_fn:} a function to apply the estimated RR transformation to new data. Takes arguments
#' \itemize{
#' \item {Y:} the data,
#' \item {z:} the z-score cutoff (defaults to 4),
#' \item {q:} the winsorizing quantile cutoff (defaults to 0.001),
#' \item {lambda:} the transformation parameter to use (defaults to the estimated one),
#' \item {T:} the transformation function family (defaults to the optimal estimated family),
#' \item {mu:} the mean to be used in the robust z-score step (re-estimates if NULL)
#' \item {sigma:} the s.d. to be used in the robust z-score step (re-estimates if NULL)
#'}
#' \item{T:} the optimal family
#' \item{T_deriv:} the derivative of the optimal family
#' \item{T_name:} name of the optimal family
#' \item{alg_control:} the parameters passed to the algorithm
#' }
#' @examples
#' Y <- rlnorm(10)%*%t(rlnorm(10))
#' rr.out <- rrscale(Y,run_parallel=FALSE)
#' Yt <- rr.out$RR
#' @importFrom foreach %do% %dopar% foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom DEoptim DEoptim DEoptim.control
rrscale <- function(Y, trans_list=list(box_cox_negative=box_cox_negative,
                                       asinh=asinh),
                    lims_list=list(box_cox_negative = c(-100,100),
                                   asinh=list(0,100)),
                  opt_control = NULL, ncores = NULL, z=4, q=0.001,
                  run_parallel = TRUE, verbose = TRUE, zeros = FALSE, opts=FALSE, seed = NULL) {
    
    set.seed(seed)

    if (is.na(zeros) | !is.logical(zeros)){
        Y[which(Y == 0)] <- zeros
    }
    
    # Find optimial parameter for each transformation family considered in trans_list
    sdns <- lapply(seq_along(trans_list), function(i) rrscale_trans(Y, trans_list[[i]], 
        lims_list[[i]], opt_control = opt_control, ncores = ncores, run_parallel = run_parallel, 
        verbose = verbose, zeros = zeros))
    names(sdns) <- names(trans_list)
    
    # Find the optimal transformation across families
    OBJ <- sapply(sdns, "[[", "objs")
    if (is.null(dim(OBJ))) 
        OBJ <- array(OBJ, c(1, length(sdns)))
    min_OBJ <- factor(apply(OBJ, 1, which.min))
    min_OBJ <- ordered(min_OBJ, levels = 1:ncol(OBJ))
    opt_trans <- which.max(table(min_OBJ))
    
    # Column-wise estimated paramters
    lambdas <- sdns[[opt_trans]]$lambdas
    
    # The median of column-wise estimates
    lambda_hat <- stats::median(lambdas, na.rm = TRUE)
    
    # Optimal transformation family
    trans <- trans_list[[opt_trans]]
    T <- trans$T
    T_opt <- T
    T_deriv <- trans$T_deriv

    # Form a function for RR that takes a matrix and returns any application of G, Z or O
    RR_fn = function(Y,G=TRUE,Z=TRUE,O=TRUE,
                     lambda=lambda_hat,T=T_opt,mu=NULL,sigma=NULL,q=0.001,z=4){

        G_fn = function(Y,lambda=lambda_hat, T=T){
            return(T(Y,lambda))
        }

        Z_fn = function(Y,mu=NULL,sigma=NULL,q=0.001){
            Y_tmp <- as.matrix(Y)
            Yw <- winsor(Y_tmp, q)
            if(is.null(mu))
                mu <- mean(Yw, na.rm = TRUE)
            Ywc <- Yw - mu
            if(is.null(sigma))
                sigma <- sqrt(mean(Ywc^2, na.rm = TRUE))
            Z <- (Y_tmp - mu)/sigma
            return(Z)
        }

        O_fn = function(Y,z=4,mu=NULL,sigma=NULL,q=0.001){
            Z = Z_fn(Y,mu=mu,sigma=sigma,q=q)
            O = Y
            O[abs(Z)>z] <- NA
            return(O)
        }

        if(G)
            Y <- G_fn(Y,lambda=lambda,T=T)
        if(Z)
            Y <- Z_fn(Y,mu=mu,sigma=sigma,q=q)
        if(O)
            Y <- O_fn(Y,z=z,mu=mu,sigma=sigma,q=q)
        return(as.matrix(Y))
    }

    RRsub = subClosure(RR_fn,c("lambda_hat","T","T_opt","q","z"))

    opt_return = NULL
    if(opts){
        opt_return = sdns
    }

    ret_list = list(opts = opt_return, # optimization output
                    pars = lambdas, # estimated columnwise parameters
                    par_hat = lambda_hat, # overall estimated parameter
                    NT = as.matrix(Y), # No transformation
                    RR = RR_fn(Y,z=z,q=q), # three-step RR transformation
                    G = RR_fn(Y,Z=FALSE,O=FALSE,z=z,q=q), # G only
                    Z = RR_fn(Y,G=FALSE,O=FALSE,z=z,q=q), # Z only
                    O = RR_fn(Y,G=FALSE,Z=FALSE,z=z,q=q), # O only 
                    rr_fn = RRsub, # rr_fn to apply to new datasets
                    T = T, # Trans family applied
                    T_deriv = T_deriv, # deriv of trans family applied
                    T_name = names(sdns)[opt_trans], # name of trans family
                    alg_control = list(trans_list, # args passed to rrscale
                                       lims_list,
                                       opt_control = opt_control,
                                       z=z,
                                       q=q,
                                       ncores = ncores,
                                       run_parallel = run_parallel, 
                                       verbose = verbose,
                                       zeros = zeros,
                                       seed = seed
                                       )
                    )
    
    return(ret_list)
}

rrscale_trans <- function(Y, tns, lims, opt_control = NULL, ncores = NULL, run_parallel = TRUE, 
    verbose = TRUE, zeros = FALSE) {
    T <- tns$T
    T_deriv <- tns$T_deriv

    if(is.logical(zeros)){
        if ((length(which(Y == 0)) > 0) && !zeros) 
            return(list(opts = NULL, lambdas = NA, objs = NA))
    }
    
    opts <- rrscale_opt(Y, lims, T, T_deriv, opt_control = opt_control, ncores = ncores, 
        run_parallel = run_parallel, verbose = verbose, zeros = zeros)
    lambdas <- sapply(opts, "[[", "est")
    objs <- sapply(opts, "[[", "obj")
    return(list(opts = opts, lambdas = lambdas, objs = objs))
}


rrscale_opt <- function(Y, lims, T, T_deriv, opt_control = NULL, ncores = NULL, run_parallel = TRUE, 
    verbose = TRUE, zeros = FALSE) {
    
    if (is.null(opt_control)) 
        opt_control <- DEoptim.control(trace = verbose, reltol = 1e-10, itermax = 10000, 
            steptol = 100)
    
    obj <- function(l, y) {
        y <- y[!is.na(y)]
        yt <- T(y, l)
        yt_deriv <- T_deriv(y, l)
        if (any(!is.finite(yt))) 
            return(Inf)
        if (diff(range(yt)) < 1e-05) 
            return(Inf)
        lsd <- log(stats::sd(yt))
        mld <- mean(log(abs(yt_deriv)))
        if (!is.finite(lsd) || !is.finite(mld)) 
            return(Inf)
        val <- lsd - mld
        return(val)
    }
    
    est_col <- function(y) {
        y <- y[!is.na(y)]
        obj_fn <- function(l) obj(l, y)
        tryCatch({
            opt <- DEoptim(fn = obj_fn, lower = lims[[1]], upper = lims[[2]], control = opt_control)
        }, error = function(e) {
            print(e)
        })
        if (!is.finite(opt$optim$bestmem)) 
            list(est = NA, obj = NA, l_tried = NA, obj_tried = NA, opt = NA)
        
        l_best <- opt$member$bestmemit
        v_best <- opt$member$bestvalit
        l_ord <- order(l_best)
        l_best <- l_best[l_ord]
        v_best <- v_best[l_ord]
        
        return(list(est = opt$optim$bestmem, obj = opt$optim$bestval, l_tried = l_best, 
            obj_tried = v_best, opt = opt))
    }
    
    operator <- `%do%`
    if (run_parallel) {
        operator <- `%dopar%`
        if (is.null(ncores)) 
            ncores <- detectCores() - 2
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
    }
    
    if (verbose) 
        dir.create(".rrscale", showWarnings = FALSE, recursive = TRUE)

    i <- 0 
    opts <- operator(foreach(i = 1:ncol(Y)), {
        fn <- NULL
        if (verbose) 
            fn <- paste0(".rrscale/", i, ".log")
        utils::capture.output({
            out <- est_col(Y[, i, drop = FALSE])
        }, file = fn)
        return(out)
    })
    
    if (run_parallel) 
        stopCluster(cl)
    
    return(opts)
}
