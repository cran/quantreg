"boot.rq"<-
function (x, y, tau = 0.5, R = 200, bsmethod = "xy", mofn = length(y), 
	  cluster = NULL, U = NULL, ...)
{
    n <- length(y)
    if(class(x) != "matrix.csr") x <- as.matrix(x)
    p <- ncol(x) 
    B <- matrix(0, R, p)
    if(tau <= 0 || tau >= 1) stop("tau outside (0,1) not allowed")
    if(length(cluster)) bsmethod <- "cluster"
    if (bsmethod == "xy") {
	if(!length(U)){
	    if(mofn < p || mofn > n) stop("mofn is out of range")
	    U <- matrix(sample(n, mofn * R, replace = TRUE), mofn, R)
	}
        B <- sqrt(mofn/n)*boot.rq.xy(x, y, U, tau)
    }   
    else if (bsmethod == "wxy") {
        if(!length(U)) U <- matrix(rexp(n * R,1), n, R)
        B <- boot.rq.wxy(x, y, U, tau)
    }   
    else if (bsmethod == "cluster"){ # Hagemann wild gradient bootstrap
	if(!length(U)){
	    if(length(cluster) != n) stop("cluster is wrong length")
	    u <- unique(cluster)
	    m <- length(u)
	    h <- c(-(sqrt(5) - 1)/2, (sqrt(5) + 1)/2)
	    hp <- c((sqrt(5) + 1)/sqrt(20), (sqrt(5) - 1)/sqrt(20))
	    U <- matrix(sample(h, R * m, prob = hp, replace = TRUE), m, R)
	    U <- apply(U, 2, function(v) v[match(cluster,u)])
	}
	r <- c(rq.fit(x, y, tau)$resid)
	psi <- (r < 0) - tau
	W <- as.matrix(t(x) %*% (U * psi))
	if(class(x) == "matrix.csr") 
	    B <- boot.rq.spwy(W, x, y, tau)
	else
	    B <- boot.rq.pwy(W, x, y, tau)
    }
    else if (bsmethod == "jack"){ # Portnoy proposal
	if(!length(U)){
	    if(missing(mofn)) mofn <- n - ceiling(2*sqrt(n))
	    if(mofn < p || mofn > n) stop("mofn is out of range")
	    U <- matrix(0, mofn, R)
	    U <- apply(U, 2, function(x) sample(n, mofn))
	}
	B <- (mofn - p + 1)/sqrt(n * (n - mofn)) * boot.rq.xy(x, y, U, tau)
    }
    else if (bsmethod == "pwy") {
	if(!length(U)){
	    U <- matrix(runif(n * R), n, R)
	}
	W <- t(x) %*% ((U < tau) - tau)
        B <- boot.rq.pwy(W, x, y, tau)
    }
    else if (bsmethod == "mcmb") {
        B <- boot.rq.mcmb(x, y, tau = tau, R = R)
    }
    else if (bsmethod == "wild") {
        n <- length(y)
        fit <- rq.fit(x, y, tau = tau)
        S <- sample(c(-2*tau,2*(1-tau)),prob = c(tau,1-tau),size = n * R, replace = TRUE)
        W <- matrix(S,n,R)
        r <- c(fit$resid)
        f0 <- akj(r,z=0)$dens 
        r <- r + hat(x) * (tau - I(r < 0))/f0 
        Y <- c(fitted(fit)) + W * abs(r)
        B <- rqs.fit(x,Y,tau = tau)
    }
    else stop("your chosen bootstrap method is not allowed")
    #cat(paste("Bootstrap standard errors based on ",R," replications"))
    list(B = B, U = U)
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
		flag = integer(R),
		coef = double(p * R),
		resid = double(m),
		integer(m),
		double((m + 5) * (p + 2)),
		double(m),
		xx = double(m * p),
		yy = double(m),
		as.integer(s),
		PACKAGE = "quantreg")
	if(sum(z$flag)>0){
		if(any(z$flag)==2) 
			warning(paste(sum(z$flag==2),"out of",R, 
				"BS replications have near singular design"))
		if(any(z$flag)==1) 
			warning(paste(sum(z$flag==1),"out of",R,"may be nonuniqu")) 
		}
	return(t(matrix(z$coef, p, R)))
}
"boot.rq.wxy"<-
function(x, y, w, tau = 0.5, tol = 0.0001)
{
#function to compute weighted bootstrap a la Bose for regression quantiles 
        x <- as.matrix(x)
        p <- ncol(x)
        n <- nrow(x)
        R <- ncol(w)
        m <- nrow(w)
        z <- .Fortran("wxy",
                as.integer(n),
                as.integer(p),
                as.integer(R),
                as.integer(m + 5),
                as.integer(p + 2),
                as.double(x),
                as.double(y),
                as.double(tau),
                as.double(tol),
                flag = integer(R),
                coef = double(p * R),
                resid = double(m),
                integer(m),
                double((m + 5) * (p + 2)),
                double(m),
                xx = double(m * p),
                yy = double(m),
                as.double(w),
                PACKAGE = "quantreg")
        if(sum(z$flag)>0){
                if(any(z$flag)==2)
                        warning(paste(sum(z$flag==2),"out of",R,
                                "BS replications have near singular design"))
                if(any(z$flag)==1)
                        warning(paste(sum(z$flag==1),"out of",R,"may be nonunique"))
                }
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
	Y <- c(y,length(y)*max(abs(y)))
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
		PACKAGE = "quantreg")
	return(t(matrix(z$coef, p, R)))
}

#################################################################
#	NB:  In an ideal world the loop would be in fortran
#################################################################


boot.rq.spwy <- function(W, a, y, tau=.5, control)
{
	y <- c(y, length(y) * max(abs(y)))
	n <- length(y)
	m <- a@dimension[2]
	a <- rbind(a, as.matrix.csr(t(rep(1,m))))
	nra <- length(a@ra)
	W <- W/tau
	k <- ncol(W) # Number of resampling realizations
	if(m != nrow(W)) 
	    stop("W row dimension not compatible with design matrix")
	if(n != a@dimension[1])
	     stop("Dimensions of design matrix and the response vector not compatible")
	for(i in 1:k){
	    a@ra[(nra - m + 1):nra] <- W[,i]
	    W[,i] <- rq.fit.sfn(a,y,tau = tau)$coef
	}
	B <- t(W)
	B
}

