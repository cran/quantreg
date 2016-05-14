#################################################################
# Interface for a sparse implementation for pwy bootstrapping
#	a = structure of the design matrix A' stored in csr format
#	y = pseudo response vector
#	W = matrix of simulated gradient vectors 
#	tau = the desired regression quantile
#	control = as for rq.fit.sfn
#################################################################


boot.rq.spwy <- function(W, a, y, tau=.5, control)
{
	Y <- -c(y, length(y) * max(abs(y)))
	n <- length(y)
	m <- a@dimension[2]
	a <- rbind(a, as(rep(1,m), "matrix.csr"))
	Wu <- W/tau
	k <- ncol(W) # Number of resampling realizations
	if(m != nrow(W)) 
	    stop("W row dimension not compatible with design matrix")
	if(n != a@dimension[1])
	     stop("Dimensions of design matrix and the response vector not compatible")
	u <- rep(1,length=n)
	x <- rep((1-tau),length=n)
	rhs = (1-tau)*c(t(a) %*% rep(1,length(y)))
	nnzdmax <- nnza <- a@ia[n+1]-1
	iwmax <- 7*m+3
	ao <- t(a)
	e <- ao %*% a
	nnzemax <- e@ia[m+1]-1
        ctrl <- sfn.control()
        if (!missing(control)) {
             control <- as.list(control)
             ctrl[names(control)] <- control
             }
	nsubmax <- ctrl$nsubmax
	tmpmax <- ctrl$tmpmax
	nnzlmax <- ctrl$nnzlmax
        if (is.null(ctrl$nsubmax)) nsubmax <- nnzemax
        if (is.null(ctrl$tmpmax)) tmpmax <- 6 * m
        if (is.null(ctrl$nnzlmax)) nnzlmax <- 4 * nnzdmax
	wwm <- vector("numeric",3*m)
	s <- u - x
	b1 <- solve(e, ao %*% y, tmpmax=tmpmax,nnzlmax=nnzlmax,nsubmax=nsubmax)
	r <- y - a %*% b1
	z <- ifelse(abs(r)<ctrl$small,(r*(r>0)+ctrl$small),r*(r>0))
	w <- z - r
	wwn <- matrix(0,n,14)
	wwn[,1] <- r
	wwn[,2] <- z
	wwn[,3] <- w
	fit <- .Fortran("spwy",
		n = as.integer(n),
		m = as.integer(m),
		k = as.integer(k),
		nnza = as.integer(nnza),
		a = as.double(a@ra),
		ja = as.integer(a@ja),
		ia = as.integer(a@ia),
		ao = as.double(ao@ra),
		jao = as.integer(ao@ja),
		iao = as.integer(ao@ia),
		nnzdmax = as.integer(nnzdmax),
		d = double(nnzdmax),
		jd = integer(nnzdmax),
		id = integer(m+1),
		dsub = double(nnzemax+1),
		jdsub = integer(nnzemax+1),
		nnzemax = as.integer(nnzemax),
		e = as.double(e@ra),
		je = as.integer(e@ja),
		ie = as.integer(e@ia),
		nsubmax = as.integer(nsubmax),
		lindx = integer(nsubmax),
		xlindx = integer(m+1),
		nnzlmax = as.integer(nnzlmax),
		lnz = double(nnzlmax),
		xlnz = integer(m+1),
		iw = integer(m*5),
		iwmax = as.integer(iwmax),
		iwork = integer(iwmax),
		xsuper = integer(m+1),
		tmpmax = as.integer(tmpmax),
		tmpvec = double(tmpmax),
		wwm = as.double(wwm),
		wwn = as.double(wwn),
		cachsz = as.integer(ctrl$cachsz),
		level = as.integer( 8 ),
		x = as.double(x),
		s = as.double(s),
		u = as.double(u),
		c = as.double(y),
		sol = as.double(Wu),
		rhs = as.double(rhs),
		small = as.double(ctrl$small),
		ierr = integer(1),
		maxiter = as.integer(ctrl$maxiter),
		time = double(7),
		PACKAGE = "quantreg")[c("sol","ierr","maxiter")]
        ierr <- fit$ierr
	if(!(ierr==0) && ctrl$warn.mesg)
            warning(sfnMessage(ierr))
	coefficients <- -fit$sol
	residuals <- -y - a %*% coefficients
	list(coefficients = coefficients,
	     residuals = residuals,
             ierr = ierr,
             it = fit$maxiter)

}
