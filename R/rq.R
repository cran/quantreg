# This are based on the rq.s file, with changes to make it work for R. 
# We try to mark the changes using some construct involving a call to 
# is.R()
# Other changes: F expanded to FALSE and T expanded to TRUE. 
#                single replaced by double. 
#                Constant PI expanded to double precision value.
#                rq: default for dual changed to TRUE (to have it coincide with doc.)
# Seems to be an error in trq, example in help won't run---neither under S-Plus. 




"qrq"<-
function(s, a)
{
#computes linearized quantiles from rq data structure
#s is the rq structure e.g. rq(x,y)
#a is a vector of quantiles required
	if(min(a) < 0 || max(a) > 1) stop("alphas out of range [0,1]")
	r <- s$sol[1,  ]
	q <- s$sol[2,  ]
	q <- c(q[1], q)
	J <- length(r)
	r <- c(0, (r[1:J - 1] + r[2:J])/2, 1)
	u <- rep(0, length(a))
	for(k in 1:length(a)) {
		i <- sum(r < a[k])
		w <- (a[k] - r[i])/(r[i + 1] - r[i])
		u[k] <- w * q[i + 1] + (1 - w) * q[i]
	}
	u
}

# rq is the function which directly calls the fortran code. 

"rq"<-
function(x, y, tau = -1, alpha = 0.10000000000000001, dual = TRUE, int = TRUE, tol
	 = 0.0001, ci = TRUE, method="score", interpolate = TRUE, tcrit = TRUE, hs=TRUE)
{
#function to compute regression quantiles
	if(!is.loaded(symbol.For("rq")))
		dyn.load("rq.dll")  #Changed(Kjetil). May need rewording for non-windows systems.
	if((tau<=0 || tau>=1) && ci){
		warning("cannot compute confidence intervals for all tau")
		ci <- FALSE
		}
	if(ci && method != "score" && method != "sparsity"){
		stop("method has to be ``score'' or ``sparsity''")
		}
	if(missing(x)) {
		x <- matrix(rep(1, length(y)), length(y))
		int <- FALSE
		if(ci && method=="score"){
			warning("method has been set to ``sparsity'' to compute confidence interval for sample quantile")
			method <- "sparsity"
			}
		}
	else{
		x <- as.matrix(x)
		if(int) 
			x <- cbind(1, x)
		else if(ncol(x)==1 && ci && method=="score"){
			warning("method has been set to ``sparsity'' to compute confidence interval when x has only one column with no intercept")
			method <- "sparsity"
			}
		}
	big <-  if( is.R() ) .Machine$double.xmax else .Machine$single.xmax
	x <- as.matrix(x)
	p <- ncol(x)
	n <- nrow(x)
	if(n != length(y))
		stop("x and y don't match number of obs.")
	nsol <- 2
	ndsol <- 2
	if(tcrit)
		cutoff <- qt(1 - alpha/2, n - p)
        else 
		cutoff <- qnorm(1 - alpha/2)
	if(ci){
		xxinv <- solve(crossprod(x))
		qn <- 1/diag(xxinv)
		}
	else
		qn <- rep(0,p)
	t.orig <- tau
	if(tau < 0 || tau > 1 || (ci && method=="sparsity")) {
		nsol <- 3 * n
		lci1 <- FALSE
		ndsol <- nsol
		tau <- -1
		}
	else {
		if(!ci)
			lci1 <- FALSE
		else lci1 <- TRUE
		}
	z <- .Fortran("rq",
		as.integer(n),
		as.integer(p),
		as.integer(n + 5),
		as.integer(p + 2),
		as.integer(p + 4),
		as.double(x),
		as.double(y),
		as.double(tau),
		as.double(tol),
		flag = as.integer(1),
		coef = double(p),
		resid = double(n),
		integer(n),
		double((n + 5) * (p + 4)),
		double(n),
		as.integer(nsol),
		as.integer(ndsol),
		sol = double((p + 2) * nsol),
		dsol = double(n * ndsol),
		lsol = as.integer(0),
		h = integer(p * nsol),
		qn = as.double(qn),
		cutoff = as.double(cutoff),
		ci = double(4 * p),
		tnmat = double(4 * p),
		as.double(big),
		as.logical(lci1))
	tau <- t.orig
	if(z$flag != 0)
		warning(switch(z$flag,"Solution may be nonunique",
		"Premature end - possible conditioning problem in x."))
	if(tau < 0 || tau > 1) {
		z$sol <- matrix(z$sol, p + 2, z$lsol)
		if(length(dimnames(x)[[2]]) == 0) {
			if(int)
				xn <- c("Intercept", paste("X", 1:(p - 1), sep
				   = ""))
			else xn <- paste("X", 1:p, sep = "")
			}
		else {
			xn <- dimnames(x)[[2]]
			if(int)
				xn[1] <- "Intercept"
			}
		dimnames(z$sol) <- list(c("Probility", "Quantile", xn), NULL)
		z$h <- matrix(z$h, p, z$lsol)
		if(dual) {
			z$dsol <- matrix(z$dsol, n, z$lsol)
			z[c("sol", "dsol", "h")]
			}
		else z[c("sol", "h")]
		}
	else {
		if(ci){
			if(method=="score") {
				if(interpolate){
					Tn <- matrix(z$tnmat, nrow = 4)
					Tci <- matrix(z$ci, nrow = 4)
					Tci[3,] <- Tci[3,] + (cutoff - abs(Tn[3,
					]))/abs(Tn[4,] - Tn[3,]) * abs(Tci[4,] - 					Tci[3,])
					Tci[2,] <- Tci[2,] - (cutoff - abs(Tn[2,
					]))/abs(Tn[1,] - Tn[2,]) * abs(Tci[1,] - 					Tci[2,])
					Tci[2,  ][is.na(Tci[2,  ])] <-  - big
					Tci[3,  ][is.na(Tci[3,  ])] <- big
					dimnames(Tci)<-list(c("Lower Bound",
					"Lower Bound","Upper Bound",
					"Upper Bound"),NULL)
					if(dual)
						return(coef = z$coef, resid = 
						z$resid, dual = z$dsol[1:n], 
						h = z$h[1:p], ci = Tci[2:3,  ])
					else 
						return(coef = z$coef, resid = 
						z$resid, h = z$h[1:p], ci = 
						Tci[2:3,  ])
					}
				else {
					Tci <- matrix(z$ci,nrow=4)
					dimnames(Tci)<-list(c("Lower Bound",
					"Lower Bound","Upper Bound",
					"Upper Bound"),NULL)
					if(dual)
						return(coef = z$coef, resid = 
						z$resid, dual = z$dsol[1:n], 
						h = z$h[1:p], ci = Tci, Tn = 
						matrix(z$tnmat, nrow = 4))
					else return(coef = z$coef, resid = 
						z$resid, h = z$h[1:p], ci = 
						Tci, Tn = matrix(z$tnmat, 
						nrow = 4))
					}
				}
			else if(method=="sparsity"){
				z$sol<-matrix(z$sol,p+2,z$lsol)
				z$dsol<-matrix(z$dsol,n,z$lsol)
				z$h <- matrix(z$h,p,z$lsol)
				se<-diag(rq.omega(z,tau,hs=hs)*xxinv)^.5
				Tci<-matrix(0,2,p)
				dimnames(Tci)<-list(c("Lower Bound",
				"Upper Bound"),NULL)
				isol <- sum(z$sol[1,]<=tau)
				Tcoef<-z$sol[3:(p+2),isol]
				Tci[1,]<-Tcoef-cutoff*se
				Tci[2,]<-Tcoef+cutoff*se
				if(dual)
					return(coef=Tcoef, resid=as.double(
					y-x%*%Tcoef), dual = z$dsol[,isol], 
					h = z$h[,isol], ci=Tci)
				else 
					return(coef=Tcoef, resid=as.double(
					y-x%*%Tcoef), h = z$h[,isol], ci=Tci)
				}
			}
		else {
			if(dual)
				return(coef = z$coef, resid = z$resid, dual = 
				z$dsol[1:n], h = z$h[1:p])
			else 
				return(coef = z$coef, resid = z$resid, h = 
				z$h[1:p])
			}
		}
}

"rq.omega"<-
function(z, tau = 0.5, hs = FALSE, tol = 0.0001)
{
#computes omega^2=tau(1-tau)/f^2(xi_a) for se of rq's
#uses linearized eqf from the full rq-process see function qrq()
#bandwidth is either Bofinger (hs=FALSE) or Hall-Sheather (hs=TRUE)
#z is the rq structure 
	dn <- dn(tau, nrow(z$dsol), hs = hs)
	low <- tau-dn
	up <- tau+dn
	if(low<0 || up>1){
		warning("bandwidth is truncated to [0,1]")
		low <- z$sol[1,2]
		up <- z$sol[1,ncol(z$sol)-1]
		}
	q <- qrq(z, c(low, up))
	shat <- diff(q)/(2 * dn)
	tau * (1 - tau) * shat^2
}

"dn"<-
function(p, n, hs = FALSE)
{
	PI <- 3.14159265358979310
	x0 <- qnorm(p)
	f0 <- (1/sqrt(2 * PI)) * exp(.Uminus((x0^2/2)))
	if(hs == TRUE)
		n^(-1/3) * qnorm(1 - 0.025000000000000001)^(2/3) * ((1.5 * f0^2
			)/(2 * x0^2 + 1))^(1/3)
	else n^-0.20000000000000001 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^
			0.20000000000000001
}

"trq"<-
function(x, y, a1 = 0.1, a2, z, int = TRUE, method = "primal", tol = 0.0001)
{
#compute trimmed regression quantiles
#z is the rq strcture 
	if(missing(a2)) a2 <- a1
	if(a1 < 0 || a2 < 0)
		stop("trimming proportion negative.")
	if(a1 + a2 - 1 > tol)
		stop("trimming proportion greater than 1.")
	if(method!="primal" && method!="dual")
		stop("invalid method: should be 'primal' or 'dual'.")
	x <- as.matrix(x)
	if(missing(z))
		z <- rq(x, y,  , int, tol)
	p <- z$sol[1,  ]
	q <- matrix(z$sol[ - c(1, 2),  ], nrow(z$sol) - 2, ncol(z$sol))
	n <- nrow(z$dsol)
	s <- NULL
	if(length(dimnames(x)[[2]]) == 0)
		dimnames(x) <- list(NULL, paste("X", 1:(nrow(q) - 1), sep = "")
			)
	if(int) {
		x <- cbind(1, x)
		dimnames(x)[[2]][1] <- "Intercept"
	}
	xbar <- apply(x, 2, "mean")
	xxinv <- solve(t(x) %*% x)
	if(abs(a1 + a2 - 1) <= tol) {
#single quantile case
		i <- sum(p < a1)
		s$coef <- q[, i]
		names(s$coef) <- dimnames(x)[[2]]
		s$resid <- y - x %*% s$coef
		PI <- 3.14159265358979310
		x0 <- qnorm(a1)
		d0 <- (1/sqrt(2 * PI)) * exp( - (x0^2/2))
		d0 <- ((4.5 * d0^4)/(2 * x0^2 + 1)^2)^0.2
		d <- d0 * (length(s$resid) - length(s$coef))^(-0.2)
		if(d > min(a1, 1 - a1))
			d <- min(a1, 1 - a1)
		s$d <- d
		i <- sum(p < a1 + d)
		j <- sum(p < a1 - d)
		shat <- as.numeric(xbar %*% t(q[, i] - q[, j]))/(2 * d)
		s$int <- int
		s$v <- a1 * (1 - a1) * shat^2
		s$cov <- s$v * xxinv
	}
	else {
#real trimming
		p1 <- p[-1]
		f <- 1/(1 - a1 - a2)
		d <- pmax((pmin(p1, 1 - a2) - c(a1, pmax(p1[ - length(p1)],
			a1))), 0)
		if(method == "primal") {
			s$coef <- q[, 1:length(p1)] %*% t(d) * f
			s$resid <- y - x %*% s$coef
			s$int <- int
		}
		else {
#Jureckova-Gutenbrunner trimmed least squares
			i <- max(1, sum(p < a1))
			g <- (z$dsol[, i + 1] - z$dsol[, i])/(p[i + 1] - p[
				i])
			wa <- z$dsol[, i] + (a1 - p[i]) * g
			j <- sum(p < 1 - a2)
			g <- (z$dsol[, j + 1] - z$dsol[, j])/(p[j + 1] - p[
				j])
			wb <- z$dsol[, j] + (1 - a2 - p[j]) * g
			wt <- wa - wb
			if(min(wt) < 0)
				warning("some weights negative!")
			s <- lsfit(x, y, abs(wt), int = FALSE)
		}
#now compute covariance matrix estimate
		mu <- xbar %*% s$coef
		v <- d %*% t((z$sol[2, 1:length(d)] - mu)^2)
		k <- qrq(z, c(a1, a2)) - mu
		v <- v + a1 * k[1]^2 + a2 * k[2]^2 + (a1 * k[1] + a2 * k[2])^2
		names(s$coef) <- dimnames(x)[[2]]
		s$v <- as.vector(f^2 * v)
		s$cov <- s$v * xxinv
	}
	s
}

"trq.print"<-
function(trq.out, digits = 4)
{
	n <- length(trq.out$resid)
	p <- length(trq.out$coef)
	options(warn = -1)
	if(trq.out$int) {
		df.num <- p - 1
		fstat <- c(trq.out$coef[-1] %*% solve(trq.out$cov[-1, -1]) %*%
			t(trq.out$coef[-1]))/df.num
	}
	else {
		df.num <- p
		fstat <- trq.out$coef %*% solve(trq.out$cov) %*% t(trq.out$coef
			)/df.num
	}
	pvalue <- 1 - pf(fstat, df.num, (n - p))
	cat(paste("Winsorized Standard Error of Regression= ", format(round(
		sqrt(trq.out$v), digits)), "\n", "N = ", format(n), 
		",  F-statistic = ", format(round(fstat, digits)), " on ",
		format(df.num), " and ", format((n - p)), " df, ", "p-value = ",
		format(round(pvalue, digits)), "\n\n", sep = ""))
	regstat <- c(sqrt(trq.out$v), n, fstat, df.num, (n - p), pvalue)
	names(regstat) <- c("rse", "n", "F.stat", "df.num", "df.den", "p.value"
		)
	err <- sqrt(diag(trq.out$cov))
	tstat <- c(trq.out$coef/err)
	tabcoef <- cbind(trq.out$coef, err, tstat, 2 * (1 - pt(abs(tstat),
		n - p)))
	dimnames(tabcoef) <- list(names(trq.out$coef), c("coef", "std.err",
		"t.stat", "p.value"))
	options(warn = 0)
	print(round(tabcoef, digits))
	invisible(list(summary = regstat, coef.table = tabcoef))
}

"ranks"<-
function(v, score = "wilcoxon")
{
	A2 <- 1
	if(score == "wilcoxon") {
		J <- ncol(v$sol)
		dt <- v$sol[1, 2:J] - v$sol[1, 1:(J - 1)]
		ranks <- as.vector((0.5 * (v$dsol[, 2:J] + v$dsol[, 1:(J - 1)]) %*%
			dt) - 0.5)
		return(ranks, A2 = 1/12)
	}
	else if(score == "normal") {
		J <- ncol(v$sol)
		dt <- v$sol[1, 2:J] - v$sol[1, 1:(J - 1)]
		dphi <- c(0, phi(qnorm(v$sol[1, 2:(J - 1)])), 0)
		dphi <- diff(dphi)
		ranks <- as.vector((((v$dsol[, 2:J] - v$dsol[, 1:(J - 1)]))) %*% (
			dphi/dt))
		return(ranks, A2)
	}
	else if(score == "sign") {
		j.5 <- sum(v$sol[1,  ] < 0.5)
		w <- (0.5 - v$sol[1, j.5])/(v$sol[1, j.5 + 1] - v$sol[1, j.5])
		r <- w * v$dsol[, j.5 + 1] + (1 - w) * v$dsol[, j.5]
		return(ranks = 2 * r - 1, A2)
	}
	else stop("invalid score function")
}

"rrs.test"<-
function(x0, x1, y, v, score = "wilcoxon")
{
	if(missing(v) || is.null(v$dsol))
		v <- rq(x0, y, dual = TRUE)
	r <- ranks(v, score)
	x1hat <- as.matrix(qr.resid(qr(cbind(1, x0)), x1))
	sn <- as.matrix(t(x1hat) %*% r$ranks)
	sn <- t(sn) %*% solve(crossprod(x1hat)) %*% sn/r$A2
	return(sn, rank = r$ranks)
}

# EXPERIMENTAL:

# Function added by Kjetil Halvorsen.
# I have followed the lm implemetation closely. 

#  "rq"<-
#  function(x, y, tau = -1, alpha = 0.10000000000000001, dual = TRUE, int = TRUE, tol
#	 = 0.0001, ci = TRUE, method="score", interpolate = TRUE, tcrit = TRUE, hs=TRUE)

"rq.formula" <- 
function(formula, data=list(), subset, na.action, tau=-1, 
         alpha = 0.10000000000000001, dual = TRUE, 
         tol = 0.0001, ci = TRUE, method="score", interpolate = TRUE, 
         tcrit = TRUE, hs=TRUE) {
# args and defaults same as in rq, + extras. Except int, which isn't necessary.
mt <- terms(formula, data=data)
mf <- mf1 <- match.call() 
mf$drop.unused.levels <- TRUE
mf$tau <- mf$alpha <- mf$dual  <- mf$tol <- mf$ci <- mf$method <- 
   mf$interpolate <- mf$tcrit <- mf$hs <- NULL
mf[[1]] <- as.name("model.frame")
mf <- eval(mf, sys.frame(sys.parent()))
xvars <- as.character(attr(mt, "variables"))[-1]
if (yvar <- attr(mt, "response") > 0)
     xvars <- xvars[-yvar]
if (length(xvars) > 0) {
   xlev <- lapply(mf[xvars], levels)
   xlev <- xlev[!sapply(xlev, is.null)]
} else xlev <- NULL
if (!is.null(model.offset(mf)))
    stop("Offset not implemented in rq.")
y <- model.response(mf, "numeric")
# Don't understand the fuss in lm about "empty.model", 
# so doesn't mimic that ...
x <- model.matrix(mt, mf)
result <- rq(x, y,  tau=tau, 
             alpha=alpha, dual=dual, int=FALSE, tol=tol, ci=ci, method=method, 
             interpolate=interpolate, tcrit=tcrit, hs=hs)
result$call <- mf1
result
}



