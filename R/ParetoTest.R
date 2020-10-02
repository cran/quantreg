ParetoTest <- function(formula, tau = 0.1, data = NULL, flavor = "Hill", 
		       m = 2, cicov = 0.9, ...) {
    if(flavor == "Hill") 
	z <-  Hill(formula, tau, data, ...)
    else if(flavor == "Pickands") 
	z <- Pickands(formula, tau, data, m = m, ...)
    else
	stop(paste("flavor", flavor, "not (yet) implemented")) 
    summary(z, ci = cicov, ...)
}

# Pickands Estimation and Inference for the Pareto Tail Exponent

Pickands <- function(formula, tau = 0.1, data = NULL, m = 2, ...){
    y <- model.response(model.frame(formula, data))
    X <- model.matrix(formula, data)
    Pickands.fit(X, y, tau, m, ...)
}
Pickands.fit <- function(X, y, tau, m, ...){
    if(tau > 0.5){ # Right tail
	X <- -X 
	y <- -y 
	tau = 1-tau
    }
    taus = c(tau, m * tau, 2 * m * tau)
    F <- rq(y ~ X - 1, tau = taus, ...)
    Xbar = colMeans(X)
    num = c(crossprod(Xbar, (F$coef[,3] - F$coef[,2])))
    denom = c(crossprod(Xbar, (F$coef[,2] - F$coef[,1])))
    xi <- (-1/log(m)) * log(num/denom)
    z <- list(xi = xi, tau = tau, m = m, F = F)
    class(z) <- "Pickands"
    z
}

summary.Pickands <- function(object, se = "boot", B = 200,  ci = 0.9, ...){
    tau <- object$tau
    xi <- object$xi
    m <- object$m
    X <- object$F$x
    y <- object$F$y
    n <- length(y)
    b1 <- object$F$coef[,1]
    bm <- object$F$coef[,2]
    Xbar <- colMeans(X)
    Xg <- X %*% (bm - b1)/c(crossprod(Xbar, bm - b1)) #aargh humbug!
    R <- rep(NA, B)
    for(i in 1:B){ # Parametric Bootstrap Default Method
	ys <- -(rexp(n)^(-xi) - 1)/xi * Xg
	R[i] <- Pickands.fit(X, ys, tau, m = m, ...)$xi
    }
    z <- rep(NA,6)
    z[1] <- xi
    z[2] <- xi - quantile(R - xi, .5, na.rm=TRUE)
    z[3] <- xi - quantile(R - xi, (1+ci)/2, na.rm=TRUE)
    z[4] <- xi - quantile(R - xi, (1-ci)/2, na.rm=TRUE)
    z[5] <- sqrt(var(R))
    z[6] <- B - sum(is.na(R))
    names(z) <- c("Estimate", "Bias-Corrected", paste("Lower ", 50*(1-ci),"%",sep=""),
	paste("Upper ", 50*(1+ci),"%",sep=""), "StdErr", "Sample Size")
    z <- list(z = z, tau = tau)
    class(z) = "summary.Pickands"
    z
}
print.Pickands <- function(x, ...){
    print(x$xi)
}
print.summary.Pickands <- function(x, ...){
    cat(paste("Pickands Estimation and Inference: tau = ", x$tau, 
	      "Effective Bootstrap Sample Size = ",x$z[6],"\n\n"))
    print(x$z[1:5])
}

if(FALSE){
# Test problem
require(quantreg)
n = 500
x = rnorm(n)
y = x + rt(n,2)
P = Pickands(y ~ x, .9)
S = summary(P)
}
# Hill Estimation and Inference for the Pareto Tail Exponent

Hill <- function(formula, tau = 0.1, data = NULL, ...){
    y <- model.response(model.frame(formula, data))
    X <- model.matrix(formula, data)
    Hill.fit(X, y, tau)
}
Hill.fit <- function(X, y, tau, ...){
    if(tau > 0.5){ # Right tail
	X <- -X 
	y <- -y 
	tau = 1-tau
    }
    F <- rq.fit(X, y, tau = tau, ...)
    yhat <- F$fitted.values
    z <- log(ifelse(y <= yhat, abs(y/yhat), 1))
    z <- z[is.finite(z)]
    z <- list(xi = sum(z)/(length(z)*tau), tau = tau, F = F)
    class(z) <- "Hill"
    z
}

summary.Hill <- function(object, se = "boot", B = 200, m = 2, ci = 0.9, ...){
    tau <- object$tau
    xi <- object$xi
    X <- object$F$x
    y <- object$F$y
    n <- length(y)
    b0 <- object$F$coef
    bm <- rq.fit(X, y, m * tau, ...)$coef 
    Xbar <- colMeans(X)
    Xg <- X %*% (bm - b0)/c(crossprod(Xbar, bm - b0)) #aargh humbug!
    R <- rep(NA, B)
    for(i in 1:B){ # Parametric Bootstrap Default Method
	ys <- -(rexp(n)^(-xi) - 1)/xi * Xg
	R[i] <- Hill.fit(X, ys, tau, ...)$xi
    }
    z <- rep(NA,6)
    z[1] <- xi
    z[2] <- xi - quantile(R - xi, .5, na.rm=TRUE)
    z[3] <- xi - quantile(R - xi, (1+ci)/2, na.rm=TRUE)
    z[4] <- xi - quantile(R - xi, (1-ci)/2, na.rm=TRUE)
    z[5] <- sqrt(var(R))
    z[6] <- B - sum(is.na(R))
    names(z) <- c("Estimate", "Bias-Corrected", paste("Lower ", 50*(1-ci),"%",sep=""),
	paste("Upper ", 50*(1+ci),"%",sep=""), "BS_StdErr", "Sample Size")
    z <- list(z = z, tau = tau)
    class(z) = "summary.Hill"
    z
}
print.Hill <- function(x, ...){
    print(x$xi)
}
print.summary.Hill <- function(x, ...){
    cat(paste("Hill Estimation and Inference: tau = ", x$tau, 
	      "Effective Bootstrap Sample Size = ",x$z[6],"\n\n"))
    print(x$z[1:5])
}
