#' a
#' @export
box_cox <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- log(Y)
        } else {
            Y <- (Y^lambda - 1)/lambda
        }
    } else {
        if (lambda == 0) {
            Y <- exp(Y)
        } else {
            Y <- (lambda * Y + 1)^(1/lambda)
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- 1/Y
    } else {
        Y <- Y^(lambda - 1)
    }
    
    return(Y)
})

#' a
#' @export
power <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        Y <- Y^lambda
    } else {
        Y <- Y^(1/lambda)
    }
    return(Y)
}, T_deriv = function(Y, lambda) {
    Y <- lambda * Y^(lambda - 1)
    return(Y)
})

#' a
#' @export
box_cox_shift <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda[1] == 0) {
            Y <- log(lambda[2] + Y)
        } else {
            Y <- ((lambda[2] + Y)^lambda[1] - 1)/lambda[1]
        }
    } else {
        if (lambd[1] == 0) {
            Y <- exp(Y) - lambda[2]
        } else {
            Y <- (lambda[1] * Y + 1)^(1/lambda) - lambda[2]
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda[1] == 0) {
        Y <- 1/(lambda[2] + Y)
    } else {
        Y <- (lambda[2] + Y)^(lambda[1] - 1)
    }
    
    return(Y)
})

#' a
#' @export
asinh <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- Y
        } else {
            Y <- asinh(lambda * Y)/lambda
        }
    } else {
        if (lambda == 0) {
            Y <- Y
        } else {
            Y <- sinh(lambda * Y)/lambda
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        if (!is.null(dim(Y))) Y <- array(1, dim = dim(Y)) else Y <- rep(1, length(Y))
    } else {
        Y <- 1/sqrt(1 + lambda^2 * Y^2)
    }
    
    return(Y)
})

#' a
#' @export
log_box_cox <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- log(log(Y))
        } else {
            Y <- log((Y^lambda - 1)/lambda)
        }
    } else {
        if (lambda == 0) {
            Y <- exp(exp(Y))
        } else {
            Y <- (lambda * exp(Y) + 1)^(1/lambda)
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- 1/(Y * log(Y))
    } else {
        Y <- (lambda * Y^(lambda - 1))/(Y^lambda - 1)
    }
    
    return(Y)
})

#' a
#' @export
box_cox_plus1 <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- log(1 + Y)
        } else {
            Y <- ((1 + Y)^lambda - 1)/lambda
        }
    } else {
        if (lambda == 0) {
            Y <- exp(Y) - 1
        } else {
            Y <- (lambda * Y + 1)^(1/lambda) - 1
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- 1/(1 + Y)
    } else {
        Y <- (1 + Y)^(lambda - 1)
    }
    
    return(Y)
})

#' a
#' @export
box_cox_negative <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        Y <- (sign(Y) * abs(Y)^lambda - 1)/lambda
    } else {
        Y <- sign(lambda * x + 1) * abs(lambda * x + 1)
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- Inf
    } else {
        Y <- abs(Y)^(lambda - 1)
    }
    
    return(Y)
})

#' a
#' @export
box_cox_plusmin <- list(T = function(Y, lambda, inverse = FALSE) {
    eps <- 1e-05
    if (min(Y) > 0) {
        l2 <- 0
    } else if (min(Y) < 0) {
        l2 <- -min(Y)
    } else if (min(Y) == 0) {
        l2 <- eps
    }
    lambda <- c(lambda, l2)
    if (!inverse) {
        if (lambda[1] == 0) {
            Y <- log(lambda[2] + Y)
        } else {
            Y <- ((lambda[2] + Y)^lambda[1] - 1)/lambda[1]
        }
    } else {
        if (lambd[1] == 0) {
            Y <- exp(Y) - lambda[2]
        } else {
            Y <- (lambda[1] * Y + 1)^(1/lambda) - lambda[2]
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    eps <- 1e-05
    if (min(Y) > 0) {
        l2 <- 0
    } else if (min(Y) < 0) {
        l2 <- -min(Y)
    } else if (min(Y) == 0) {
        l2 <- eps
    }
    lambda <- c(lambda, l2)
    if (lambda[1] == 0) {
        Y <- 1/(lambda[2] + Y)
    } else {
        Y <- (lambda[2] + Y)^(lambda[1] - 1)
    }
    
    return(Y)
})

#' a
#' @export
box_cox_exp <- list(T = function(Y, lambda, inverse = FALSE) {
    if (!inverse) {
        if (lambda == 0) {
            Y <- Y
        } else {
            Y <- (exp(lambda * Y) - 1)/lambda
        }
    } else {
        if (lambda == 0) {
            Y <- Y
        } else {
            Y <- log(lambda * Y + 1)/lambda
        }
    }
    
    return(Y)
}, T_deriv = function(Y, lambda) {
    if (lambda == 0) {
        Y <- Y
    } else {
        Y <- exp(lambda * Y)
    }
    
    return(Y)
})

#' List possible transformations. 
#' @export
list_transformations <- function() {
    return(list(box_cox = box_cox, power = power, box_cox_shift = box_cox_shift, 
        asinh = asinh, log_box_cox = log_box_cox, box_cox_plus1 = box_cox_plus1, 
        box_cox_negative = box_cox_negative, box_cox_plusmin = box_cox_plusmin, box_cox_exp = box_cox_exp))
}
