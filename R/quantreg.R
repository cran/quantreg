".First.lib" <-
function(lib, pkg) {
   library.dynam("quantreg", pkg, lib)
   print("quantreg library loaded")}

"bandwidth.rq" <-
function(p, n, hs = TRUE, alpha = 0.05)
{
	# Bandwidth selection for sparsity estimation two flavors:
	#	Hall and Sheather(1988, JRSS(B)) rate = O(n^{-1/3})
	#	Bofinger (1975, Aus. J. Stat)  -- rate = O(n^{-1/5})
	# Generally speaking, default method, hs=TRUE is preferred.
	PI <- pi
	x0 <- qnorm(p)
	f0 <- (1/sqrt(2 * PI)) * exp( - ((x0^2/2)))
	if(hs == TRUE)
		n^(-1/3) * qnorm(1 - alpha/2)^(2/3) * ((1.5 * f0^2)/(2 * x0^
			2 + 1))^(1/3)
	else n^-0.2 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^ 0.2
}


"plot.rq.process" <-
# Function to plot estimated quantile regression  process
function(x, nrow = 3, ncol = 2, ...)
{
	tdim <- dim(x$sol)
	p <- tdim[1] - 3
	m <- tdim[2]
	par(mfrow = c(nrow, ncol))
	ylab <- dimnames(x$sol)[[1]]
	for(i in 1:p) {
		plot(x$sol[1,], x$sol[3 + i,  ], xlab = "tau", ylab = ylab[3 +
			i], type = "l")
	}
}

"print.rqs" <-
function (x, ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    coef <- coef(x)
    cat("\nCoefficients:\n")
    print(coef, ...)
    rank <- x$rank
    nobs <- nrow(residuals(x))
    p <- nrow(coef)
    rdf <- nobs - p
    cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
    if (!is.null(attr(x, "na.message"))) 
        cat(attr(x, "na.message"), "\n")
    invisible(x)
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

"print.summary.rqs" <- 
function(x, ...) lapply(x,print.summary.rq)

"print.summary.rq" <-
function(x, digits = max(5, .Options$digits - 2), ...)
{
	cat("\nCall: ")
	dput(x$call)
	coef <- x$coef
	df <- x$df
	rdf <- x$rdf
	tau <- x$tau
	cat("\ntau: ")
	print(format(round(tau,digits = digits)), quote = FALSE, ...)
	cat("\nCoefficients:\n")
	print(format(round(coef, digits = digits)), quote = FALSE, ...)
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
		return(list(ranks=ranks, A2=A2))
	}
	else if(score == "normal") {
		J <- ncol(v$sol)
		dt <- v$sol[1, 2:J] - v$sol[1, 1:(J - 1)]
		dphi <- c(0, dnorm(qnorm(v$sol[1, 2:(J - 1)])), 0)
		dphi <- diff(dphi)
		ranks <- as.vector((((v$dsol[, 2:J] - v$dsol[, 1:(J - 1)]))) %*%
			(dphi/dt))
		return(list(ranks=ranks, A2=A2))
	}
	else if(score == "sign") {
		j.5 <- sum(v$sol[1,  ] < 0.5)
		w <- (0.5 - v$sol[1, j.5])/(v$sol[1, j.5 + 1] - v$sol[1, j.5])
		r <- w * v$dsol[, j.5 + 1] + (1 - w) * v$dsol[, j.5]
		ranks <- 2 * r - 1
		return(list(ranks=ranks, A2=A2))
	}
	else if(score == "tau") {
		j.tau <- sum(v$sol[1,  ] < tau)
		w <- (tau - v$sol[1, j.tau])/(v$sol[1, j.tau + 1] - v$sol[
			1, j.tau])
		r <- w * v$dsol[, j.tau + 1] + (1 - w) * v$dsol[, j.tau]
		ranks <- r - (1-tau)
        	A2 <- tau * (1-tau)
		return(list(ranks=ranks, A2=A2))
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
		return(list(ranks=ranks, A2=A2))
	}
	else stop("invalid score function")
}

"rq" <-
function (formula, tau = 0.5, data, weights, na.action, method = "br", 
    model = TRUE, contrasts = NULL, ...) 
{
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
    mf <- mf[c(1,m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    if(method == "model.frame")return(mf)
    mt <- attr(mf, "terms")
    weights <- model.weights(mf)
    Y <- model.response(mf)
    X <- model.matrix(mt, mf, contrasts)
    if(length(tau)>1){
	if(any(tau <0) || any(tau >1)) stop("invalid taus")
	coef <- matrix(0,ncol(X),length(tau))
	resid <- matrix(0,nrow(X),length(tau))
	for(i in 1:length(tau)){
	    z <- {if (length(weights)) 
                 	rq.wfit(X, Y, tau = tau[i], weights, method, ...)
                  else rq.fit(X, Y, tau = tau[i], method, ...)
    		 }
	    coef[,i] <- z$coefficients 
	    resid[,i] <- z$residuals
	   }
	taulabs <- paste("tau=",format(round(tau,3)))
	dimnames(coef) <- list(dimnames(X)[[2]],taulabs)
	dimnames(resid) <- list(dimnames(X)[[1]],taulabs)
	fit <- list(coefficients = coef, residuals = resid)
	class(fit) <- "rqs"
	}
    else{
        process <- (tau < 0 || tau > 1)
        fit <- {
            if (length(weights)) 
                rq.wfit(X, Y, tau = tau, weights, method, ...)
            else rq.fit(X, Y, tau = tau, method, ...)
           }
        class(fit) <- ifelse(process, "rq.process", "rq")
	}
    fit$formula <- formula
    fit$terms <- mt
    fit$call <- call
    fit$tau <- tau
    attr(fit, "na.message") <- attr(m, "na.message")
    if(model) fit$model <- mf
    fit
}
"rq.fit" <-
function(x, y, tau = 0.5, method = "br", ...)
{
	fit <- switch(method,
		fn = rq.fit.fn(x, y, tau = tau, ...),
		fnb = rq.fit.fnb(x, y, tau = tau, ...),
		fnc = rq.fit.fnc(x, y, tau = tau, ...),
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

# Function to compute regression quantiles using original simplex approach
# of Barrodale-Roberts/Koenker-d'Orey.  There are several options.
# The options are somewhat different than those available for the Frisch-
# Newton version of the algorithm, reflecting the different natures of the
# problems typically solved.  Succintly BR for "small" problems, FN for
# "large" ones.  Obviously, these terms are conditioned by available hardware.
#
# Basically there are two modes of use:
# 1.  For Single Quantiles:
#
#       if tau is between 0 and 1 then only one quantile solution is computed.
#
#       if ci = FALSE  then just the point estimate and residuals are returned
#		If the column dimension of x is 1 then ci is set to FALSE since
#		since the rank inversion method has no proper null model.
#       if ci = TRUE  then there are two options for confidence intervals:
#
#               1.  if iid = TRUE we get the original version of the rank
#                       inversion intervals as in Koenker (1994)
#               2.  if iid = FALSE we get the new version of the rank inversion
#                       intervals which accounts for heterogeneity across
#                       observations in the conditional density of the response.
#                       The theory of this is described in Koenker-Machado(1999)
#               Both approaches involve solving a parametric linear programming
#               problem, the difference is only in the factor qn which
#               determines how far the PP goes.  In either case one can
#               specify two other options:
#                       1. interp = FALSE returns two intervals an upper and a
#                               lower corresponding to a level slightly
#                               above and slightly below the one specified
#                               by the parameter alpha and dictated by the
#                               essential discreteness in the test statistic.
#				interp = TRUE  returns a single interval based on
#                               linear interpolation of the two intervals
#                               returned:  c.values and p.values which give
#                               the critical values and p.values of the
#                               upper and lower intervals. Default: interp = TRUE.
#                       2.  tcrit = TRUE uses Student t critical values while
#                               tcrit = FALSE uses normal theory ones.
# 2. For Multiple Quantiles:
#
#       if tau < 0 or tau >1 then it is presumed that the user wants to find
#       all of the rq solutions in tau, and the program computes the whole
#	quantile regression solution as a process in tau, the resulting arrays
#	containing the primal and dual solutions, betahat(tau), ahat(tau)
#       are called sol and dsol.  These arrays aren't printed by the default
#       print function but they are available as attributes.
#       It should be emphasized that this form of the solution can be
#	both memory and cpu quite intensive.  On typical machines it is
#	not recommended for problems with n > 10,000.
#	In large problems a grid of solutions is probably sufficient.
#
rq.fit.br <- 
function (x, y, tau = 0.5, alpha = 0.1, ci = FALSE, iid = TRUE, interp = TRUE, tcrit = TRUE) 
{
    tol <- .Machine$double.eps^(2/3)
    eps <- tol
    big <- .Machine$double.xmax
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    nsol <- 2
    ndsol <- 2
    if (tau < 0 || tau > 1) {
        nsol <- 3 * n
        ndsol <- 3 * n
        lci1 <- FALSE
        qn <- rep(0, p)
        cutoff <- 0
        tau <- -1
    }
    else {
        if (p == 1) 
            ci <- FALSE
        if (ci) {
            lci1 <- TRUE
            if (tcrit) 
                cutoff <- qt(1 - alpha/2, n - p)
            else cutoff <- qnorm(1 - alpha/2)
            if (!iid) {
                h <- bandwidth.rq(tau, n, hs = TRUE)
                bhi <- rq.fit.br(x, y, tau + h, ci = FALSE)
                bhi <- coefficients(bhi)
                blo <- rq.fit.br(x, y, tau - h, ci = FALSE)
                blo <- coefficients(blo)
                dyhat <- x %*% (bhi - blo)
                if (any(dyhat <= 0)) {
                  pfis <- (100 * sum(dyhat <= 0))/n
                  warning(paste(pfis, "percent fis <=0"))
                }
                f <- pmax(eps, (2 * h)/(dyhat - eps))
                qn <- rep(0, p)
                for (j in 1:p) {
                  qnj <- lm(x[, j] ~ x[, -j] - 1, weights = f)$resid
                  qn[j] <- sum(qnj * qnj)
                }
            }
            else qn <- 1/diag(solve(crossprod(x)))
        }
        else {
            lci1 <- FALSE
            qn <- rep(0, p)
            cutoff <- 0
        }
    }
    z <- .Fortran("rqbr", as.integer(n), as.integer(p), as.integer(n + 
        5), as.integer(p + 3), as.integer(p + 4), as.double(x), 
        as.double(y), as.double(tau), as.double(tol), flag = as.integer(1), 
        coef = double(p), resid = double(n), integer(n), double((n + 
            5) * (p + 4)), double(n), as.integer(nsol), as.integer(ndsol), 
        sol = double((p + 3) * nsol), dsol = double(n * ndsol), 
        lsol = as.integer(0), h = integer(p * nsol), qn = as.double(qn), 
        cutoff = as.double(cutoff), ci = double(4 * p), tnmat = double(4 * 
            p), as.double(big), as.logical(lci1), PACKAGE = "quantreg")
    if (z$flag != 0) 
        warning(switch(z$flag, "Solution may be nonunique", "Premature end - possible conditioning problem in x"))
    if (tau < 0 || tau > 1) {
        sol <- matrix(z$sol[1:((p + 3) * z$lsol)], p + 3)
        dsol <- matrix(z$dsol[1:(n * z$lsol)], n)
        vnames <- dimnames(x)[[2]]
        dimnames(sol) <- list(c("tau", "Qbar", "Obj.Fun", vnames), 
            NULL)
        return(list(sol = sol, dsol = dsol))
    }
    if (!ci) {
        coef <- z$coef
        names(coef) <- dimnames(x)[[2]]
        return(list(coefficients = coef, x = x, y = y, residuals = z$resid))
    }
    if (interp) {
        Tn <- matrix(z$tnmat, nrow = 4)
        Tci <- matrix(z$ci, nrow = 4)
        Tci[3, ] <- Tci[3, ] + (abs(Tci[4, ] - Tci[3, ]) * (cutoff - 
            abs(Tn[3, ])))/abs(Tn[4, ] - Tn[3, ])
        Tci[2, ] <- Tci[2, ] - (abs(Tci[1, ] - Tci[2, ]) * (cutoff - 
            abs(Tn[2, ])))/abs(Tn[1, ] - Tn[2, ])
        Tci[2, ][is.na(Tci[2, ])] <- -big
        Tci[3, ][is.na(Tci[3, ])] <- big
        coefficients <- cbind(z$coef, t(Tci[2:3, ]))
        vnames <- dimnames(x)[[2]]
        cnames <- c("coefficients", "lower bd", "upper bd")
        dimnames(coefficients) <- list(vnames, cnames)
        residuals <- z$resid
        return(list(coefficients = coefficients, residuals = residuals))
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
        p.values <- if (tcrit) 
            matrix(pt(c.values, n - p), ncol = 4)
        else matrix(pnorm(c.values), ncol = 4)
        dimnames(p.values) <- list(vnames, cnames[-1])
        list(coefficients = coefficients, residuals = residuals, 
            c.values = c.values, p.values = p.values)
    }
}
"rq.fit.fn" <-
function (x, y, tau = 0.5, beta = 0.99995, eps = 1e-06)
{
    n <- length(y)
    p <- ncol(x)
    if (n != nrow(x))
        stop("x and y don't match n")
    if (tau < eps || tau > 1 - eps)
        stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
    rhs <- (1 - tau) * apply(x, 2, sum)
    d   <- rep(1,n)
    u   <- rep(1,n)
    wn <- rep(0,10*n)
    wn[1:n] <- (1-tau) #initial value of dual solution
    z <- .Fortran("rqfn", as.integer(n), as.integer(p), a = as.double(t(as.matrix(x))),
        c = as.double(-y), rhs = as.double(rhs), d = as.double(d),as.double(u),
        beta = as.double(beta), eps = as.double(eps), 
        wn = as.double(wn), wp = double((p + 3) * p), aa = double(p *
            p), it.count = integer(2), info = integer(1),PACKAGE= "quantreg")
    if (z$info != 0)
        stop(paste("Error info = ", z$info, "in stepy: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    residuals <- y - x %*% coefficients
    list(coefficients=coefficients, tau=tau, residuals=residuals)
}

"rq.fit.fnb" <-
function (x, y, tau = 0.5, beta = 0.99995, eps = 1e-06)
{
    n <- length(y)
    p <- ncol(x)
    if (n != nrow(x))
        stop("x and y don't match n")
    if (tau < eps || tau > 1 - eps)
        stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
    rhs <- (1 - tau) * apply(x, 2, sum)
    d   <- rep(1,n)
    u   <- rep(1,n)
    wn <- rep(0,10*n)
    wn[1:n] <- (1-tau) #initial value of dual solution
    z <- .Fortran("rqfnb", as.integer(n), as.integer(p), a = as.double(t(as.matrix(x))),
        c = as.double(-y), rhs = as.double(rhs), d = as.double(d),as.double(u),
        beta = as.double(beta), eps = as.double(eps), 
        wn = as.double(wn), wp = double((p + 3) * p), 
        it.count = integer(2), info = integer(1),PACKAGE= "quantreg")
    if (z$info != 0)
        stop(paste("Error info = ", z$info, "in stepy: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    residuals <- y - x %*% coefficients
    list(coefficients=coefficients, tau=tau, residuals=residuals)
}

"rq.fit.fnc" <-
function (x, y, R, r, tau = 0.5, beta = 0.9995, eps = 1e-06) 
{
    n1 <- length(y)
    n2 <- length(r)
    p <- ncol(x)
    if (n1 != nrow(x)) 
        stop("x and y don't match n1")
    if (n2 != nrow(R)) 
        stop("R and r don't match n2")
    if (p != ncol(R)) 
        stop("R and x don't match p")
    if (tau < eps || tau > 1 - eps) 
        stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
    rhs <- (1 - tau) * apply(x, 2, sum)
    u <- rep(1, n1) #upper bound vector
    wn1 <- rep(0, 9 * n1)
    wn1[1:n1] <- (1 - tau) #store the values of x1
    wn2 <- rep(0, 6 * n2) 
    wn2[1:n2] <- 1 #store the values of x2
    z <- .Fortran("rqfnc", as.integer(n1), as.integer(n2), as.integer(p), 
	a1 = as.double(t(as.matrix(x))), c1 = as.double(-y), 
	a2 = as.double(t(as.matrix(R))), c2 = as.double(-r), 
	rhs = as.double(rhs), d1 = double(n1), d2 = double(n2),
	as.double(u), beta = as.double(beta), eps = as.double(eps), 
	wn1 = as.double(wn1), wn2 = as.double(wn2), wp = double((p + 3) * p), 
	it.count = integer(2), info = integer(1),PACKAGE= "quantreg")
    if (z$info != 0) 
        stop(paste("Error info = ", z$info, "in stepy2: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    residuals <- y - x %*% coefficients
    it.count <- z$it.count
    list(coefficients=coefficients, tau=tau, residuals=residuals)
}


"rq.fit.pfn" <-
# This is an implementation (purely in R) of the preprocessing phase
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
function(x, y, tau = 0.5,  Mm.factor = 0.8,
	max.bad.fixup = 3, eps = 1e-6)
{
	#rq function for n large --
	n <- length(y)
	if(nrow(x) != n)
		stop("x and y don't match n")
	if(tau < 0 | tau > 1)
		stop("tau outside (0,1)")
	p <- ncol(x)
	m <- round(((p + 1) * n)^(2/3))
	not.optimal <- TRUE
	while(not.optimal) {
		if(m < n)
			s <- sample(n, m)
		else {
			z <- rq.fit.fn(x, y, tau = tau,  eps = eps)
			list(coef = z$coef)
		}
		xx <- x[s,  ]
		yy <- y[s]
		z <- rq.fit.fn(xx, yy, tau = tau,  eps = eps)
		xxinv <- solve(chol(crossprod(xx)))
		band <- sqrt(((x %*% xxinv)^2) %*% rep(1, p))
		#sqrt(h<-ii)
		r <- y - x %*% z$coef
		M <- Mm.factor * m
		lo.q <- max(1/n, tau - M/(2 * n))
		hi.q <- min(tau + M/(2 * n), (n - 1)/n)
		kappa <- quantile(r/pmax(eps, band), c(lo.q, hi.q))
		sl <- r < band * kappa[1]
		su <- r > band * kappa[2]
		bad.fixup <- 0
		while(not.optimal & (bad.fixup < max.bad.fixup)) {
			xx <- x[!su & !sl,  ]
			yy <- y[!su & !sl]
			if(any(sl)) {
				glob.x <- c(t(x[sl,  , drop = FALSE]) %*% rep(
					1, sum(sl)))
				glob.y <- sum(y[sl])
				xx <- rbind(xx, glob.x)
				yy <- c(yy, glob.y)
			}
			if(any(su)) {
				ghib.x <- c(t(x[su,  , drop = FALSE]) %*% rep(
					1, sum(su)))
				ghib.y <- sum(y[su])
				xx <- rbind(xx, ghib.x)
				yy <- c(yy, ghib.y)
			}
			z <- rq.fit.fn(xx, yy, tau = tau,  eps = eps)
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
			else not.optimal <- FALSE
		}
	}
	coefficients <- b
	names(coefficients) <- dimnames(x)[[2]]
	residuals <- y - x %*% b
	return(list(coefficients=coefficients, tau=tau, 
		residuals=residuals))
}

"rq.wfit" <-
function(x, y, tau = 0.5, weights, method = "br",  ...)
{
	if(any(weights < 0))
		stop("negative weights not allowed")
	contr <- attr(x, "contrasts")
	x <- x * weights
	y <- y * weights
	fit <- switch(method,
		fn = rq.fit.fn(x, y, tau = tau, ...),
		br = rq.fit.br(x, y, tau = tau, ...),
		fnc = rq.fit.fnc(x, y, tau = tau, ...), 
                pfn = rq.fit.pfn(x, y, tau = tau, ...), {
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
	return(list(sn=sn, ranks=ranks))
}

"summary.rqs" <-
function (object, ...) {
        taus <- object$tau
        xsum <- as.list(taus)
        for(i in 1:length(taus)){
                xi <- object
                xi$coefficients <- xi$coefficients[,i]
                xi$residuals <- xi$residuals[,i]
                xi$tau <- xi$tau[i]
                class(xi) <- "rq"
                xsum[[i]] <- summary(xi, ...)
                }
        class(xsum) <- "summary.rqs"
   	mt <- terms(object)
    	m <- model.frame(object)
    	y <- model.response(m)
    	x <- model.matrix(mt,m,contrasts = object$contrasts)
	xsum$olscoefs <- coef(summary(lm(y~x)))
        xsum
        }

"summary.rq" <-
# This function provides  methods for summarizing the output of the
# rq command. In this instance, "summarizing" means essentially provision
# of either standard errors, or confidence intervals for the rq coefficents.
# Since the preferred method for confidence intervals is currently the
# rank inversion method available directly from rq() by setting ci=TRUE, with br=TRUE.
# these summary methods are intended primarily for comparison purposes
# and for use on large problems where the parametric programming methods
# of rank inversion are prohibitively memory/time consuming.  Eventually
# iterative versions of rank inversion should be developed that would
# employ the Frisch-Newton approach.
#
# Object is the result of a call to rq(), and the function returns a
# table of coefficients, standard errors, "t-statistics", and p-values, and, if
# covariance=TRUE a structure describing the covariance matrix of the coefficients,
# i.e. the components of the Huber sandwich.
#
# There are five options for "se":
#
#	1.  "rank" strictly speaking this doesn't produce a "standard error"
#		at all instead it produces a coefficient table with confidence
#		intervals for the coefficients based on inversion of the
#		rank test described in GJKP and Koenker (1994).  
#	2.  "iid" which presumes that the errors are iid and computes
#		an estimate of the asymptotic covariance matrix as in KB(1978).
#	3.  "nid" which presumes local (in tau) linearity (in x) of the
#		the conditional quantile functions and computes a Huber
#		sandwich estimate using a local estimate of the sparsity.
#	4.  "ker" which uses a kernel estimate of the sandwich as proposed
#		by Powell.
#	5.  "boot" which uses a bootstrap method:
#		"xy"	uses xy-pair method 
#		"pwy"	uses the parzen-wei-ying method 
#		"mcmb"	uses the Markov chain marginal bootstrap method 
#
#
function (object, se = "nid", covariance = TRUE, hs = TRUE, ...) 
{
    mt <- terms(object)
    m <- model.frame(object)
    y <- model.response(m)
    x <- model.matrix(mt,m,contrasts = object$contrasts)
    wt <- model.weights(m)
    tau <- object$tau
    eps <- .Machine$double.eps^(2/3)
    coef <- coefficients(object)
    if (is.matrix(coef)) 
        coef <- coef[, 1]
    vnames <- dimnames(x)[[2]]
    resid <- object$residuals
    n <- length(resid)
    p <- length(coef)
    rdf <- n - p
    if (!is.null(wt)) {
        resid <- resid * wt
        x <- x * wt
        y <- y * wt
    }
    if (missing(se)) {
        if (n < 1001) 
            se <- "rank"
        else se <- "nid"
    }
    if (se == "rank") {
        f <- rq.fit.br(x, y, tau = tau, ci = TRUE, ...)
    }
    if (se == "iid") {
        xxinv <- diag(p)
        xxinv <- backsolve(qr(x)$qr[1:p, 1:p], xxinv)
        xxinv <- xxinv %*% t(xxinv)
        pz <- sum(abs(resid) < eps)
        h <- max(p + 1, ceiling(n * bandwidth.rq(tau, n, hs = hs)))
        ir <- (pz + 1):(h + pz + 1)
        ord.resid <- sort(resid[order(abs(resid))][ir])
        xt <- ir/(n - p)
        sparsity <- rq(ord.resid ~ xt)$coef[2, 1]
        cov <- sparsity^2 * xxinv * tau * (1 - tau)
        serr <- sqrt(diag(cov))
    }
    else if (se == "nid") {
        h <- bandwidth.rq(tau, n, hs = hs)
        if (tau + h > 1) 
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0) 
            stop("tau - h < 0:  error in summary.rq")
        bhi <- rq.fit.fn(x, y, tau = tau + h)$coef
        blo <- rq.fit.fn(x, y, tau = tau - h)$coef
        dyhat <- x %*% (bhi - blo)
        if (any(dyhat <= 0)) 
            warning(paste(sum(dyhat <= 0), "non-positive fis"))
        f <- pmax(0, (2 * h)/(dyhat - eps))
        fxxinv <- diag(p)
        fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p], fxxinv)
        fxxinv <- fxxinv %*% t(fxxinv)
        cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*% 
            fxxinv
        serr <- sqrt(diag(cov))
    }
    else if (se == "ker") {
        h <- bandwidth.rq(tau, n, hs = hs)
        if (tau + h > 1) 
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0) 
            stop("tau - h < 0:  error in summary.rq")
        uhat <- c(y - x %*% coef)
        h <- (qnorm(tau + h) - qnorm(tau - h))*
		min(sqrt(var(uhat)), ( quantile(uhat,.75)- quantile(uhat, .25))/1.34 )
        f <- dnorm(uhat/h)/h
        fxxinv <- diag(p)
        fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p], fxxinv)
        fxxinv <- fxxinv %*% t(fxxinv)
        cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*% 
            fxxinv
        serr <- sqrt(diag(cov))
    }
    else if (se == "boot") {
        B <- boot.rq(x, y, tau, ...)
        cov <- cov(B)
        serr <- sqrt(diag(cov))
    }
    if( se == "rank"){
	coef <- f$coef
	}
    else {
    	coef <- array(coef, c(p, 4))
    	dimnames(coef) <- list(vnames, c("Value", "Std. Error", "t value", 
             "Pr(>|t|)"))
    	coef[, 2] <- serr
    	coef[, 3] <- coef[, 1]/coef[, 2]
    	coef[, 4] <- if (rdf > 0) 
			2 * (1 - pt(abs(coef[, 3]), rdf))
    		     else NA
	}
    object <- object[c("call", "terms")]
    if (covariance == TRUE) {
        object$cov <- cov
        if (se %in% c("nid", "ker")) {
            object$Hinv <- fxxinv
            object$J <- crossprod(x)
        }
        else if (se == "boot") {
            object$B <- B
        }
    }
    object$coefficients <- coef
    object$rdf <- rdf
    object$tau <- tau
    class(object) <- "summary.rq"
    object
}

"akj" <-
function(x, z = seq(min(x), max(x),  , 2 * length(x)), 
	p = rep(1, length(x))/ length(x), h = 0, alpha = 0.5, kappa = 0.9, 
	iker1 = 1, iker2 = 1)
{
        nx <- length(x)
        nz <- length(z)
        x <- sort(x)
        A <- .Fortran("akj",
                as.double(x),
                as.double(z),
                as.double(p),
                as.integer(iker1),
                dens = double(nz),
                psi = double(nz),
                score = double(nz),
                as.integer(nx),
                as.integer(nz),
                h = as.double(h),
                as.double(alpha),
                as.double(kappa),
                double(nx),
		PACKAGE = "quantreg")
        dens <- A$dens
        psi <- A$psi
        score <- A$score
        h <- A$h
        return(list(dens=dens, psi=psi, score=score, h=h))
}
"lm.fit.recursive" <-
function(X, y, int = TRUE)
{
	if(int)
		X <- cbind(1, X)
	p <- ncol(X)
	n <- nrow(X)
	D <- qr(X[1:p,  ])
	A <- qr.coef(D, diag(p))
	A[is.na(A)] <- 0
	A <- crossprod(t(A))
	Ax <- rep(0, p)
	b <- matrix(0, p, n)
	b[, p] <- qr.coef(D, y[1:p])
	b[is.na(b)] <- 0
	z <- .Fortran( "rls",
		as.integer(n),
		as.integer(p),
		as.double(t(X)),
		as.double(y),
		b = as.double(b),
		as.double(A),
		as.double(Ax), 
		PACKAGE = "quantreg")
	bhat <- matrix(z$b, p, n)
	return(bhat)
}
