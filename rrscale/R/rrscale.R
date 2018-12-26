#' a
#' @export
#' @param Y Data matrix to be transformed.
#' @param trans_list List of transformations to be considered. See function list_transformations. 
#' @param lims_list List of optimization limits for each transformation from trans_list. 
#' @param opt_control Optional optimization controlling parameters for DEoptim control argument. 
#' @param ncores Number of cores to use if running parallel. 
#' @param run_parallel TRUE/FALSE run column-wise in parallel. 
#' @param verbose TRUE/FALSE save optimization output in local directory '.rrscale'
#' @param zeros How to deal with zeros in the data set. If set to FALSE the algorithm will fail if it encounters a zero. If set to 'adjust' then the zeros are replaced by 1E-10. If set to 'remove' then the zeros are set to NA. 
#' @param seed Sets the seed before running any other analyses. 
#' @return A list of output:
#' opts: the optimization output for all transformation families and all columns
#' pars: the optimal parameters for each column for the optimal family
#' par_hat: the estimated optimal paramter
#' NT: the original data
#' RR: the robust-rescaled data
#' G: gaussianized data
#' Z: robust z-transformed data
#' O: data with outliers removed
#' transform: a function that takes the data and performs the RR transformation given the estimated values of the transformation parameter, mean, and sd of the data.
#' T: the optimal family
#' T_deriv: the derivative of the optimal family
#' T_name: name of the optimal family
#' alg_control: the parameters passed to the algorithm
#' @importFrom foreach %do% %dopar% foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom DEoptim DEoptim DEoptim.control
rrscale <- function(Y, trans_list, lims_list, opt_control = NULL, ncores = NULL, 
    run_parallel = TRUE, verbose = TRUE, zeros = FALSE, seed = NULL) {
    
    if (!is.null(seed)) {
        set.seed(seed)
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
    T_deriv <- trans$T_deriv
    
    Yt_global <- T(Y, lambda_hat)
    Yt_global <- as.matrix(Yt_global)
    
    Yw_global <- winsor(Yt_global, 0.001)
    mu_global <- mean(Yw_global, na.rm = TRUE)
    Ywc_global <- Yw_global - mu_global
    norm_global <- sqrt(mean(Ywc_global^2, na.rm = TRUE))
    
    transform <- function(Y, center = TRUE, scale = TRUE, rm_q = 4, re_estimate = FALSE) {
        
        Yt <- T(Y, lambda_hat)
        Yt <- as.matrix(Yt)
        
        if (re_estimate) {
            Yw <- winsor(Yt, 0.001)
            mu <- mean(Yw, na.rm = TRUE)
            Ywc <- Yw - mu
            norm <- sqrt(mean(Ywc^2, na.rm = TRUE))
        } else {
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
    
    Yt <- transform(Y)
    
    subTransform <- subClosure(transform, c("mu_global", "norm_global", "winsor"))
    
    Y_tmp <- as.matrix(Y)
    Yw <- winsor(Y_tmp, 0.001)
    mu <- mean(Yw, na.rm = TRUE)
    Ywc <- Yw - mu
    norm <- sqrt(mean(Ywc^2, na.rm = TRUE))
    Y_tmp <- (Y_tmp - mu)/norm
    Z <- Y_tmp
    
    O <- Y
    O[abs(Z) > 4] <- NA
    
    return(list(opts = sdns, pars = lambdas, par_hat = lambda_hat, NT = as.matrix(Y), 
        RR = as.matrix(Yt), G = as.matrix(transform(Y, center = FALSE, scale = FALSE, 
            rm_q = Inf)), Z = as.matrix(Z), O = as.matrix(O), transform = subTransform, 
        T = T, T_deriv = T_deriv, T_name = names(sdns)[opt_trans], alg_control = list(trans_list, 
            lims_list, opt_control = opt_control, ncores = ncores, run_parallel = run_parallel, 
            verbose = verbose, zeros = zeros, seed = seed)))
}

rrscale_trans <- function(Y, tns, lims, opt_control = NULL, ncores = NULL, run_parallel = TRUE, 
    verbose = TRUE, zeros = FALSE) {
    T <- tns$T
    T_deriv <- tns$T_deriv
    
    if ((length(which(Y == 0)) > 0) && !zeros) 
        return(list(opts = NULL, lambdas = NA, objs = NA))
    
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
    
    if (zeros == "adjust") 
        Y[which(Y == 0)] <- 1e-10
    if (zeros == "remove") {
        Y[which(Y == 0)] <- NA
    }
    
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
