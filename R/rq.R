"Tn.test" <-
function(x0, x1, y, m = 100)
{
	fit <- rq(cbind(x0, x1), y, 0.5, ci = T)
	p <- ncol(cbind(1, x0, x1))
	cilen <- fit$ci[2, p] - fit$ci[1, p]
	v <- fit$coef[p] + (0.5 * cilen * ( - m:m))/m
	Tn <- v
	for(i in 1:length(v)) {
		Tn[i] <- rrs.test(x0, x1, y - x1 * v[i], score = "sign")$sn
	}
	return(v, Tn, ci = fit$ci)
}

"bandwidth.rq" <-
function(p, n, hs = T, alpha = 0.05)
{
	# Bandwidth selection for sparsity estimation two flavors:
	#	Hall and Sheather(1988, JRSS(B)) rate = O(n^{-1/3})
	#	Bofinger (1975, Aus. J. Stat)  -- rate = O(n^{-1/5})
	# Generally speaking, default method, hs=T is preferred.
	PI <- pi
	x0 <- qnorm(p)
	f0 <- (1/sqrt(2 * PI)) * exp( - ((x0^2/2)))
	if(hs == T)
		n^(-1/3) * qnorm(1 - alpha/2)^(2/3) * ((1.5 * f0^2)/(2 * x0^
			2 + 1))^(1/3)
	else n^-0.2 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^ 0.2
}


"plot.rq.process" <-
# Function to plot estimated quantile regression  process
function(object, nrow = 3, ncol = 2)
{
	tdim <- dim(object$sol)
	p <- tdim[1] - 3
	m <- tdim[2]
	par(mfrow = c(nrow, ncol))
	x <- object$sol[1,  ]
	ylab <- dimnames(object$sol)[[1]]
	for(i in 1:p) {
		plot(x, object$sol[3 + i,  ], xlab = "tau", ylab = ylab[3 +
			i], type = "l")
	}
}

"plot.table.rq" <-
# Function to plot estimated quantile regression parameters in table.rq
function(object, nrow = 3, ncol = 2)
{
	tdim <- dim(object$a)
	p <- tdim[1]
	m <- tdim[2]
	k <- tdim[3]
	par(mfrow = c(nrow, ncol))
	x <- object$taus
	for(i in 1:p) {
		ylab <- dimnames(object$a)[[1]]
		plot(rep(x, 2), object$a[i,  , 2:3], xlab = "tau", ylab = ylab[
			i], type = "n")
		points(x, object$a[i,  , 1], pch = "o")
		lines(x, object$a[i,  , 1])
		lines(x, object$a[i,  , 2], lty = 2)
		lines(x, object$a[i,  , 3], lty = 2)
	}
}

"print.rq" <-
function(x, ...)
{
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}
	coef <- coef(x)
	cat("\nCoefficients:\n")
	print(coef, ...)
	rank <- x$rank
	nobs <- length(residuals(x))
	if(is.matrix(coef))
		p <- dim(coef)[1]
	else p <- length(coef)
	rdf <- nobs - p
	cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
	if(!is.null(attr(x, "na.message")))
		cat(attr(x, "na.message"), "\n")
	invisible(x)
}

"print.summary.rq" <-
function(x, digits = max(5, .Options$digits - 2), ...)
{
	cat("\nCall: ")
	dput(x$call)
	coef <- x$coef
	df <- x$df
	rdf <- x$rdf
	cat("\nCoefficients:\n")
	print(format(round(coef, digits = digits)), quote = F, ...)
	invisible(x)
}

"ranks" <-
function(v, score = "wilcoxon", tau = 0.5)
{
	A2 <- 1
	if(score == "wilcoxon") {
		J <- ncol(v$sol)
		dt <- v$sol[1, 2:J] - v$sol[1, 1:(J - 1)]
		ranks <- as.vector((0.5 * (v$dsol[, 2:J] + v$dsol[, 1:(J - 1)]) %*%
			dt) - 0.5)
		A2 <- 1/12
		return(ranks, A2)
	}
	else if(score == "normal") {
		J <- ncol(v$sol)
		dt <- v$sol[1, 2:J] - v$sol[1, 1:(J - 1)]
		dphi <- c(0, phi(qnorm(v$sol[1, 2:(J - 1)])), 0)
		dphi <- diff(dphi)
		ranks <- as.vector((((v$dsol[, 2:J] - v$dsol[, 1:(J - 1)]))) %*%
			(dphi/dt))
		return(ranks, A2)
	}
	else if(score == "sign") {
		j.5 <- sum(v$sol[1,  ] < 0.5)
		w <- (0.5 - v$sol[1, j.5])/(v$sol[1, j.5 + 1] - v$sol[1, j.5])
		r <- w * v$dsol[, j.5 + 1] + (1 - w) * v$dsol[, j.5]
		ranks <- 2 * r - 1
		return(ranks, A2)
	}
	else if(score == "tau") {
		j.tau <- sum(v$sol[1,  ] < tau)
		w <- (tau - v$sol[1, j.tau])/(v$sol[1, j.tau + 1] - v$sol[
			1, j.tau])
		r <- w * v$dsol[, j.tau + 1] + (1 - w) * v$dsol[, j.tau]
		ranks <- 2 * r - 1
		return(ranks, A2)
	}
	else if(score == "interquartile") {
		j.25 <- sum(v$sol[1,  ] < 0.25)
		w <- (0.25 - v$sol[1, j.25])/(v$sol[1, j.25 + 1] - v$sol[1,
			j.25])
		r.25 <- w * v$dsol[, j.25 + 1] + (1 - w) * v$dsol[, j.25]
		j.75 <- sum(v$sol[1,  ] < 0.75)
		w <- (0.75 - v$sol[1, j.75])/(v$sol[1, j.75 + 1] - v$sol[1,
			j.75])
		r.75 <- w * v$dsol[, j.75 + 1] + (1 - w) * v$dsol[, j.75]
		ranks <- 0.5 + r.75 - r.25
		A2 <- 1/4
		return(ranks, A2)
	}
	else stop("invalid score function")
}

"rq" <-
function(formula, tau = 0.5, data, weights, na.action, method = "br", contrasts
	 = NULL, ...)
{
	call <- match.call()
	m <- match.call(expand = F)
	m$method <- m$model <- m$x <- m$y <- m$contrasts <- m$tau <- m$... <-
		NULL
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.frame(sys.parent()))
	if(method == "model.frame")
		return(m)
	Terms <- attr(m, "terms")
	weights <- model.extract(m, weights)
	Y <- model.extract(m, response)
	X <- model.matrix(Terms, m, contrasts)
	process <- (tau < 0 || tau > 1)
	fit <- {
		if(length(weights))
			rq.wfit(X, Y, tau = tau, weights, method, ...)
		else rq.fit(X, Y, tau = tau, method, ...)
	}
	fit$terms <- Terms
	fit$call <- call
	fit$tau <- tau
	attr(fit, "na.message") <- attr(m, "na.message")
	class(fit) <- ifelse(process, "rq.process", "rq")
	fit
}

"rq.fit" <-
function(x, y, tau = 0.5, method = "br", ...)
{
	#if(!is.numeric(x)) stop("model matrix must be numeric")
	#if(!is.numeric(y)) stop("response must be numeric")
	fit <- switch(method,
		fn = rq.fit.fn(x, y, tau = tau, ...),
		pfn = rq.fit.pfn(x, y, tau = tau, ...),
		br = rq.fit.br(x, y, tau = tau, ...),
		{
			what <- paste("rq.fit.", method, sep = "")
			if(exists(what, mode = "function"))
				(get(what, mode = "function"))(x, y, ...)
			else stop(paste("unimplemented method:", method))
		}
		)
	fit$contrasts <- attr(x, "contrasts")
	fit
}

"rq.fit.br" <-
# Function to compute regression quantiles using original simplex approach
# of Barrodale-Roberts/Koenker-d'Orey.  There are several options.
# The options are somewhat different than those available for the Frisch-
# Newton version of the algorithm, reflecting the different natures of the
# problems typically solved.  Succintly BR for "small" problems, FN for
# "large" ones.  Obviously, these terms are conditioned by available hardware.
#
# Basically there are two modes of use:
# For Single Quantiles:
#
#       if tau is between 0 and 1 then only one quantile solution is computed.
#
#       if ci = F  then just the point estimate and residuals are returned
#       if ci = T  then there are two options for confidence intervals:
#
#               1.  if iid = T we get the original version of the rank
#                       inversion intervals as in Koenker (1994)
#               2.  if iid = F we get the new version of the rank inversion
#                       intervals which accounts for heterogeneity across
#                       observations in the conditional density of the response.#                       The theory of this is described in Koenker-Machado(1999)#
#               Both approaches involve solving a parametric linear programming
#               problem, the difference is only in the factor qn which
#               determines how far the PP goes.  In either case one can
#               specify two other options:
#                       1. interp = F returns two intervals an upper and a
#                               lower corresponding to a level slightly
#                               above and slightly below the one specified
#                               by the parameter alpha and dictated by the
#                               essential discreteness in the test statistic.
#				interp = T  returns a single interval based on
#                               linear interpolation of the two intervals
#                               returned:  c.values and p.values which give
#                               the critical values and p.values of the
#                               upper and lower intervals. Default: interp = T.
#                       2.  tcrit = T uses Student t critical values while
#                               tcrit = F uses normal theory ones.
# 2. For Multiple Quantiles:
#
#       if tau < 0 or tau >1 then it is presumed that the user wants to find
#       all of the rq solutions in tau, and the program computes the whole
#	quantile regression solution as a process in tau, the resulting arrays
#	containing the primal and dual solutions, betahat(tau), ahat(tau)
#       are called sol and dsol.  These arrays aren't printed by the default
#       function print.rq but they are available as attributes.
#       It should be emphasized that this form of the solution can be
#	both memory and cpu quite intensive.  On typical machines it is
#	not recommended for problems with n > 10,000.
#	In large problems a grid of solutions is probably sufficient.
#
function(x, y, tau = 0.5, alpha = 0.10000000000000001, ci = T, iid = T, interp
	 = T, tcrit = T)
{
	tol <- .Machine$double.eps^(2/3)
	eps <- tol
	big <- .Machine$double.xmax
	x <- as.matrix(x)
	p <- ncol(x)
	n <- nrow(x)
	nsol <- 2
	ndsol <- 2
	if(tau < 0 || tau > 1) {
		nsol <- 3 * n
		ndsol <- 3 * n
		lci1 <- F
		qn <- rep(0, p)
		cutoff <- 0
		tau <- -1
	}
	else {
		if(ci) {
			lci1 <- T
			if(tcrit)
				cutoff <- qt(1 - alpha/2, n - p)
			else cutoff <- qnorm(1 - alpha/2)
			if(!iid) {
				h <- bandwidth.rq(tau, n, hs = T)
				bhi <- rq.fit.br(x, y, tau + h, ci = F)
				bhi <- coefficients(bhi)
				blo <- rq.fit.br(x, y, tau - h, ci = F)
				blo <- coefficients(blo)
				dyhat <- x %*% (bhi - blo)
				if(any(dyhat <= 0)) {
					pfis <- (100 * sum(dyhat <= 0))/n
					warning(paste(pfis, "percent fis <=0"))
				}
				f <- pmax(eps, (2 * h)/(dyhat - eps))
				qn <- rep(0, p)
				for(j in 1:p) {
					qnj <- lm(x[, j] ~ x[,  - j] - 1,
						weights = f)$resid
					qn[j] <- sum(qnj * qnj)
				}
			}
			else qn <- 1/diag(solve(crossprod(x)))
		}
		else {
			lci1 <- F
			qn <- rep(0, p)
			cutoff <- 0
		}
	}
	assign <- attr(x, "assign")
	z <- .Fortran("rqbr",
		as.integer(n),
		as.integer(p),
		as.integer(n + 5),
		as.integer(p + 3),
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
		sol = double((p + 3) * nsol),
		dsol = double(n * ndsol),
		lsol = as.integer(0),
		h = integer(p * nsol),
		qn = as.double(qn),
		cutoff = as.double(cutoff),
		ci = double(4 * p),
		tnmat = double(4 * p),
		as.double(big),
		as.logical(lci1))
	if(z$flag != 0)
		warning(switch(z$flag,
			"Solution may be nonunique",
			"Premature end - possible conditioning problem in x"))
	if(tau < 0 || tau > 1) {
		sol <- matrix(z$sol[1:((p + 3) * z$lsol)], p + 3)
		dsol <- matrix(z$dsol[1:(n * z$lsol)], n)
		vnames <- dimnames(x)[[2]]
		dimnames(sol) <- list(c("tau", "Qbar", "Obj.Fun", vnames),
			NULL)
		return(sol, dsol)
	}
	if(!ci) {
		coefficients <- z$coef
		residuals <- z$resid
		return(coefficients, x, y, residuals)
	}
	if(interp) {
		Tn <- matrix(z$tnmat, nrow = 4)
		Tci <- matrix(z$ci, nrow = 4)
		Tci[3,  ] <- Tci[3,  ] + (abs(Tci[4,  ] - Tci[3,  ]) * (cutoff -
			abs(Tn[3,  ])))/abs(Tn[4,  ] - Tn[3,  ])
		Tci[2,  ] <- Tci[2,  ] - (abs(Tci[1,  ] - Tci[2,  ]) * (cutoff -
			abs(Tn[2,  ])))/abs(Tn[1,  ] - Tn[2,  ])
		Tci[2,  ][is.na(Tci[2,  ])] <-  - big
		Tci[3,  ][is.na(Tci[3,  ])] <- big
		coefficients <- cbind(z$coef, t(Tci[2:3,  ]))
		vnames <- dimnames(x)[[2]]
		cnames <- c("coefficients", "lower bd", "upper bd")
		dimnames(coefficients) <- list(vnames, cnames)
		residuals <- z$resid
		return(coefficients, x, y, residuals, assign)
	}
	else {
		Tci <- matrix(z$ci, nrow = 4)
		coefficients <- cbind(z$coef, t(Tci))
		residuals <- z$resid
		vnames <- dimnames(x)[[2]]
		cnames <- c("coefficients", "lower bound", "Lower Bound",
			"upper bd", "Upper Bound")
		dimnames(coefficients) <- list(vnames, cnames)
		c.values <- t(matrix(z$tnmat, nrow = 4))
		c.values <- c.values[, 4:1]
		dimnames(c.values) <- list(vnames, cnames[-1])
		p.values <- matrix(pt(c.values, n - p), ncol = 4)
		dimnames(p.values) <- list(vnames, cnames[-1])
		return(coefficients, x, y, residuals, c.values, p.values)
	}
}

"rq.fit.fn" <-
function(x, y, tau = 0.5, int = F, beta = 0.99995, eps = 1e-6)
{
	#rq function using Frisch-Newton interior point algorithm
	n <- length(y)
	if(int)
		x <- cbind(1, x)
	p <- ncol(x)
	if(n != nrow(x))
		stop("x and y don't match n")
	if(tau < eps || tau > 1 - eps)
		stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
	z <- .Fortran("rqfn",
		as.integer(n),
		as.integer(p),
		a = as.double(t(as.matrix(x))),
		c = as.double(y),
		rhs = double(p),
		d = double(n),
		beta = as.double(beta),
		eps = as.double(eps),
		tau = as.double(tau),
		wn = double(10 * n),
		wp = double((p + 3) * p),
		aa = double(p * p),
		it.count = integer(2),
		info = integer(1))
	if(z$info != 0)
		stop(paste("Error", z$info, "in stepy: singular design"))
	coefficients <- z$wp[1:p]
	names(coefficients) <- dimnames(x)[[2]]
	residuals <- y - x %*% coefficients
	return(coefficients, x, y, tau, residuals)
}

"rq.fit.pfn" <-
# This is an implementation (purely in S) of the preprocessing phase
# of the rq algorithm described in Portnoy and Koenker, Statistical
# Science, (1997) 279-300.  In this implementation it can be used
# as an alternative method for rq() by specifying method="pfn"
# It should probably be used only on very large problems and then
# only with some caution.  Very large in this context means roughly
# n > 100,000.  The options are described in the paper, and more
# explicitly in the code.  Again, it would be nice perhaps to have
# this recoded in a lower level language, but in fact this doesn't
# seem to make a huge difference in this case since most of the work
# is already done in the rq.fit.fn calls.
#
function(x, y, tau = 0.5, int = F, Mm.factor = 0.8,
	max.bad.fixup = 3, eps = 1e-6)
{
	#rq function for n large --
	n <- length(y)
	if(int)
		x <- cbind(1, x)
	if(nrow(x) != n)
		stop("x and y don't match n")
	if(tau < 0 | tau > 1)
		stop("tau outside (0,1)")
	p <- ncol(x)
	m <- round(((p + 1) * n)^(2/3))
	not.optimal <- T
	while(not.optimal) {
		if(m < n)
			s <- sample(n, m)
		else {
			z <- rq.fit.fn(x, y, tau = tau, int = F)
			return(coef = z$coef)
		}
		xx <- x[s,  ]
		yy <- y[s]
		z <- rq.fit.fn(xx, yy, tau = tau, int = F)
		xxinv <- solve(chol(crossprod(xx)))
		band <- sqrt(((x %*% xxinv)^2) %*% rep(1, p))
		#sqrt(h_ii)
		r <- y - x %*% z$coef
		M <- Mm.factor * m
		lo.q <- max(1/n, tau - M/(2 * n))
		hi.q <- min(tau + M/(2 * n), (n - 1)/n)
		kappa <- kuantile(r/pmax(eps, band), c(lo.q, hi.q))
		sl <- r < band * kappa[1]
		su <- r > band * kappa[2]
		bad.fixup <- 0
		while(not.optimal & (bad.fixup < max.bad.fixup)) {
			xx <- x[!su & !sl,  ]
			yy <- y[!su & !sl]
			if(any(sl)) {
				glob.x <- c(t(x[sl,  , drop = F]) %*% rep(
					1, sum(sl)))
				glob.y <- sum(y[sl])
				xx <- rbind(xx, glob.x)
				yy <- c(yy, glob.y)
			}
			if(any(su)) {
				ghib.x <- c(t(x[su,  , drop = F]) %*% rep(
					1, sum(su)))
				ghib.y <- sum(y[su])
				xx <- rbind(xx, ghib.x)
				yy <- c(yy, ghib.y)
			}
			z <- rq.fit.fn(xx, yy, tau = tau, int = F)
			b <- z$coef
			r <- y - x %*% b
			su.bad <- (r < 0) & su
			sl.bad <- (r > 0) & sl
			if(any(c(su.bad, sl.bad))) {
				if(sum(su.bad | sl.bad) > 0.10000000000000001 *
					M) {
					warning("Too many fixups:  doubling m")
					m <- 2 * m
					break
				}
				su <- su & !su.bad
				sl <- sl & !sl.bad
				bad.fixup <- bad.fixup + 1
			}
			else not.optimal <- F
		}
	}
	coefficients <- b
	names(coefficients) <- dimnames(x)[[2]]
	residuals <- y - x %*% b
	return(coefficients, x, y, tau, residuals)
}

"rq.wfit" <-
function(x, y, tau = 0.5, weights, method, int = T, ...)
{
	if(!is.numeric(x))
		stop("model matrix must be numeric")
	if(!is.numeric(y))
		stop("response must be numeric")
	if(any(weights < 0))
		stop("negative weights not allowed")
	contr <- attr(x, "contrasts")
	x <- x * weights
	y <- y * weights
	fit <- switch(method,
		fn = rq.fit.fn(x, y, tau = tau, ...),
		br = rq.fit.br(x, y, tau = tau, ...),
		{
			what <- paste("rq.fit.", method, sep = "")
			if(exists(what, mode = "function"))
				(get(what, mode = "function"))(x, y, ...)
			else stop(paste("unimplemented method:", method))
		}
		)
	fit$contrasts <- attr(x, "contrasts")
	fit
}

"rrs.test" <-
function(x0, x1, y, v, score = "wilcoxon")
{
	if(missing(v) || is.null(v$dsol))
		v <- rq(y ~ x0, tau = -1)
	r <- ranks(v, score)
	x1hat <- as.matrix(qr.resid(qr(cbind(1, x0)), x1))
	sn <- as.matrix(t(x1hat) %*% r$ranks)
	sn <- t(sn) %*% solve(crossprod(x1hat)) %*% sn/r$A2
	ranks <- r$ranks
	return(sn, ranks)
}

"summary.rq" <-
# This is a preliminary method for summarizing the output of the
# rq command eventually some bootstrapping strategies should be
# added.  In this instance, "summarizing" means essentially provision
# of either standard errors, or confidence intervals for the rq coefficents.
# Since the preferred method for confidence intervals is currently the
# rank inversion method available directly from rq() by setting ci=T, with br=T.
# these summary methods are intended primarily for comparison purposes
# and for use on large problems where the parametric programming methods
# of rank inversion are prohibitively memory/time consuming.  Eventually
# iterative versions of rank inversion should be developed that would
# employ the Frisch-Newton approach.
#
# Object is the result of a call to rq(), and the function returns a
# table of coefficients, standard errors, "t-statistics", and p-values, and, if
# covariance=T a structure describing the covariance matrix of the coefficients,
# i.e. the components of the Huber sandwich.
#
# There are three options for "se":
#
#	1.  "iid" which presumes that the errors are iid and computes
#		an estimate of the asymptotic covariance matrix as in KB(1978).
#	2.  "nid" which presumes local (in tau) linearity (in x) of the
#		the conditional quantile functions and computes a Huber
#		sandwich estimate using a local estimate of the sparsity.
#	3.  "ker" which uses a kernel estimate of the sandwich as proposed
#		by Powell.
# See the inference chapter of the putative QR book for further details.
#
#
function(object, se = "nid", covariance = T)
{
	x <- object$x
	y <- object$y
	tau <- object$tau
	eps <- .Machine$double.eps^(2/3)
	wt <- object$weights
	coef <- coefficients(object)
	if(is.matrix(coef))
		coef <- coef[, 1]
	vnames <- dimnames(x)[[2]]
	resid <- object$residuals
	n <- length(resid)
	p <- length(coef)
	rdf <- n - p
	if(!is.null(wt)) {
		resid <- resid * wt
		x <- x * wt
		y <- y * wt
	}
	#quick and dirty se's in three flavors: iid, nid, and ker
	if(se == "iid") {
		xxinv <- diag(p)
		xxinv <- backsolve(qr(x)$qr[1:p, 1:p], xxinv)
		xxinv <- xxinv %*% t(xxinv)
		pz <- sum(abs(resid) < eps)
		h <- max(p + 1, ceiling(n * bandwidth.rq(tau, n, hs = T)))
		ir <- (pz + 1):(h + pz + 1)
		ord.resid <- sort(resid[order(abs(resid))][ir])
		sparsity <- l1fit(ir/(n - p), ord.resid)$coef[2]
		cov <- sparsity^2 * xxinv * tau * (1 - tau)
		serr <- sqrt(diag(cov))
	}
	else if(se == "nid") {
		h <- bandwidth.rq(tau, n, hs = T)
		bhi <- rq.fit.fn(x, y, tau = tau + h, int = F)$coef
		blo <- rq.fit.fn(x, y, tau = tau - h, int = F)$coef
		dyhat <- x %*% (bhi - blo)
		if(any(dyhat <= 0))
			warning(paste(sum(dyhat <= 0), "non-positive fis"))
		f <- pmax(0, (2 * h)/(dyhat - eps))
		fxxinv <- diag(p)
		fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p], fxxinv)
		fxxinv <- fxxinv %*% t(fxxinv)
		cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*% fxxinv
		serr <- sqrt(diag(cov))
	}
	else if(se == "ker") {
		h <- bandwidth.rq(tau, n, hs = T)
		h <- qnorm(tau + h) -  qnorm(tau - h)
		uhat <- c(y - x %*% coef)
		f <- dnorm(uhat/h)/h
		fxxinv <- diag(p)
		fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p], fxxinv)
		fxxinv <- fxxinv %*% t(fxxinv)
		cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*% fxxinv
		serr <- sqrt(diag(cov))
	}
	coef <- array(coef, c(p, 4))
	dimnames(coef) <- list(vnames, c("Value", "Std. Error", "t value",
		"Pr(>|t|)"))
	coef[, 2] <- serr
	coef[, 3] <- coef[, 1]/coef[, 2]
	coef[, 4] <- if(rdf > 0) 2 * (1 - pt(abs(coef[, 3]), rdf)) else NA
#	object <- object[c("call", "terms", "assign")]
	object <- object[c("call", "terms")]
	if(covariance == T) {
		object$cov <- cov
		if(se != "iid") {
			object$Hinv <- fxxinv
			object$J <- crossprod(x)
		}
	}
	object$coefficients <- coef
	object$rdf <- rdf
	class(object) <- "summary.rq"
	object
}

"table.rq" <-
# This is a function to produce a table of estimated coefficients
# for a quantile regression problem.  The default method
# is "br", i.e. simplex, with confidence intervals for the coefficients
# estimated by the rank inversion method.  For large problems it may be
# prudent to adopt an alternative strategy.  The function returns an
# array of dimension (p,m,3) where p is the dimension of betahat, m is
# the length of the vector of taus, and 3 accounts for coefs,upper,lower
# or coefs,se,t-stat depending on the method invoked.
#
function(formula, taus = c(0.25, 0.5, 0.75), method = "br", ...)
{
	m <- length(taus)
	tab <- NULL
	for(i in 1:m) {
		fit <- rq(formula, taus[i], method = method)
		tab <- rbind(tab, coefficients(fit))
	}
	p <- nrow(tab)/m
	a <- array(tab, dim = c(p, m, 3))
	vnames <- dimnames(coefficients(fit))[[1]]
	ctypes <- c("coefs", "upper ci limit", "lower ci limit")
	dimnames(a) <- list(vnames, paste("tau=", taus), ctypes)
	tab <- list(a = a, taus = taus)
	class(tab) <- "table.rq"
	invisible(tab)
}

"translate" <-
function(text, old, new, multichar = F)
{
	if(length(old) > 1 || (nchar(old) != nchar(new)))
		multichar <- T
	if(length(old) > 1 && (length(new) > 1 & length(new) != length(old)))
		stop("old and new must have same lengths or new must have 1 element"
			)
	if(multichar)
		command <- paste("sed", paste("-e 's/", old, "/", new, "/g'",
			sep = "", collapse = " "))
	else command <- paste("tr \"", old, "\" \"", new, "\"", sep = "")
	k <- system(command, text)
	if(is.matrix(text))
		k <- matrix(k, nrow = nrow(text))
	k
}

