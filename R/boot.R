"boot.rq"<-
function (x, y, tau = 0.5, R = 200, bsmethod = "xy", mofn = n, ...)
{
    n <- length(y)
    x <- as.matrix(x)
    p <- ncol(x) 
    B <- matrix(0, R, p)
    if (bsmethod == "xy") {
	if(mofn < p || mofn > n) stop("mofn is out of range")
        s <- matrix(sample(n, mofn * R, replace = TRUE), mofn, R)
        B <- sqrt(mofn/n)*boot.rq.xy(x, y, s, tau)
    }   
    else if (bsmethod == "pwy") {
        U <- t(x) %*% matrix(((runif(n * R) > tau) - tau), n,
            R)
        B <- boot.rq.pwy(U, x, y, tau)
    }
    else if (bsmethod == "mcmb") {
        B <- boot.rq.mcmb(x, y, tau = tau, R = R)
    }
    else stop("your chosen bootstrap method is not allowed")
    #cat(paste("Bootstrap standard errors based on ",R," replications"))
    B
}
"boot.rq.mcmb" <-
function (x, y, tau = 0.5, R = 200) 
{
	n <- length(y)
	p <- ncol(x)
	if(n < 2000)
		fit <- rq(y~x - 1, tau = tau, ci = FALSE)
	else
		fit <- rq(y~x - 1, tau = tau, method = "fn")
	coef <- fit$coef
	svdx <- svd(x)
	condx <- svdx$d[1]/svdx$d[p]
	Ainv <- svdx$v %*% diag(svdx$d) %*% t(svdx$v) 
	coefTilda <-  Ainv %*% coef
	A <- svdx$v %*% diag(1/svdx$d) %*% t(svdx$v) 
	r <- fit$resid
    	psi <- signr <- sign(r)
    	psi[signr > 0] <- tau
    	psi[signr < 0] <- tau - 1
	psimat <- matrix(psi, nrow = n, ncol = p, byrow = FALSE)
	x <- x %*% A
	ZTilda <- x * psimat
	sumxij <- apply(x, 2, sum)
	sumabsxij <- apply(abs(x), 2, sum)
	zstar <- .C("bootnp", as.double(t(x)), as.double(y), 
		as.double(tau), as.double(coefTilda), 
		as.double(t(A)), as.double(ZTilda), 
		as.double(sumxij), as.double(sumabsxij), 
		as.integer(n), as.integer(p), success = as.integer(1), 
		theta = as.double(rep(0, R * p + p), as.integer(c(p, R + 1))), 
		as.integer(R), PACKAGE = "quantreg")
	if (zstar$success == 0) 
		return(list(success = 0))
	else{
		B <- matrix(zstar$theta,p,R+1)[,-1]
		B <- t(A %*% B)
		}
	B
}

"boot.rq.xy"<-
function(x, y, s, tau = 0.5, tol = 0.0001)
{
#function to compute xypairs bootstrap for regression quantiles 
#stripped down for monte-carlo purposes
	x <- as.matrix(x)
	p <- ncol(x)
	n <- nrow(x)
	R <- ncol(s)
	m <- nrow(s)
	z <- .Fortran("xys",
		as.integer(m),
		as.integer(n),
		as.integer(p),
		as.integer(R),
		as.integer(m + 5),
		as.integer(p + 2),
		as.double(x),
		as.double(y),
		as.double(tau),
		as.double(tol),
		flag = as.integer(1),
		coef = double(p * R),
		resid = double(m),
		integer(m),
		double((m + 5) * (p + 2)),
		double(m),
		as.integer(1),
		sol = double((p + 2)),
		dsol = double(m),
		lsol = as.integer(0),
		xx = double(m * p),
		yy = double(m),
		as.integer(s),
		PACKAGE = "quantreg")
	return(t(matrix(z$coef, p, R)))
}
"boot.rq.pwy"<-
function(U,X, y, tau = 0.5, tol=1e-4)
{
#resampling method of parzen,wei,ying for quantile regression
#NB. x should be full design matrix including intercept
	n <- length(y)
	p <- ncol(X)
	R <- ncol(U)
	Y <- c(y,500000)
	x <- rbind(X,0)
	xu <- t(U)/tau
	n <- n+1
	z<-.Fortran("pwy",
		as.integer(n),
		as.integer(p),
		as.integer(R),
		as.integer(n + 5),
		as.integer(p + 2),
		as.double(xu),
		as.double(x),
		as.double(Y),
		as.double(tau),
		as.double(tol),
		flag = as.integer(1),
		coef = double(p * R),
		resid = double(n),
		integer(n),
		double((n + 5) * (p + 2)),
		double(n),
		as.integer(1),
		sol = double((p + 2)),
		dsol = double(n),
		lsol = as.integer(0),
		PACKAGE = "quantreg")
	return(t(matrix(z$coef, p, R)))
}
