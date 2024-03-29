"boot.rq"<-
function (x, y, tau = 0.5, R = 200, bsmethod = "xy", mofn = length(y), 
	  coef = NULL, blbn = NULL, cluster = NULL, U = NULL, ...)
{
    n <- length(y)
    if(class(x)[1] != "matrix.csr") x <- as.matrix(x)
    p <- ncol(x) 
    B <- matrix(0, R, p)
    if(tau < 0 || tau > 1) stop("tau outside (0,1) not allowed")
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
    else if (bsmethod == "pxy") {
        if (!length(U)) {
            if (mofn < p || mofn > n) stop("mofn is out of range")
            U <- matrix(sample(n, mofn * R, replace = TRUE), mofn, R)
            }
        B <- sqrt(mofn/n) * boot.rq.pxy(x, y, U, tau, coef, ...)
    }
    else if (bsmethod == "pwxy") {
        B <- sqrt(mofn/n) * boot.rq.pwxy(x, y, tau, coef, ...)$B
    }
    else if (bsmethod == "BLB") {
        B <- boot.rq.pwxy(x, y, tau, coef, ...)$B
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
	if(class(x)[1] == "matrix.csr") 
	    r <- c(rq.fit.sfn(x, y, tau)$resid)
	else
	    r <- c(rq.fit.fnb(x, y, tau)$resid)
	psi <- (r < 0) - tau
	W <- as.matrix(t(x) %*% (U * psi))
	if(class(x)[1] == "matrix.csr") 
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
	W <- as.matrix(t(x) %*% ((U < tau) - tau))
	if(class(x)[1] == "matrix.csr") 
	    B <- boot.rq.spwy(W, x, y, tau)
	else
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
	if(inherits(x, "matrix.csr"))
	    H <- diag(x %*% solve(t(x) %*% x) %*% t(x))
	else
	    H <- hat(x, intercept = FALSE) 
        r <- r + H * (tau - I(r < 0))/f0 
        Y <- c(fitted(fit)) + W * abs(r)
        B <- rqs.fit(x,Y,tau = tau)
    }
    else stop("your bootstrap method is not implemented")
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
		as.integer(R))
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
		as.integer(s))
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
                as.double(w))
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
	if(inherits(X, "matrix.csr"))
	    x <- rbind(X, t(as.matrix.csr(rep(0, p))))
	else
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
		double(n))
	return(t(matrix(z$coef, p, R)))
}

#################################################################
#	NB:  In an ideal world the loop below would be in fortran, in a
#		a really ideal world all of this resampling code
#		would be able to adapt to the "fn" or "sfn" etc
#		method used by the original fitting scheme.
#		My initial attempt to allow control to be passed
#		from the initial fit of the coefficients was a failure.
#		There was a "memory not mapped" segfault presumably
#		due to the fact that the pwy scheme has one extra
#		observation.  Who would have guessed things were so
#		delicately tuned.  This should be fixed when I'm in
#		a more delicate frame of mind.  Meanwhile, on Andreas
#		test problem using sfn is substantially faster 26s
#		vs 148s for fnb, and even more for br as I used earlier.
#################################################################


boot.rq.spwy <- function(W, a, y, tau=.5)
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
# Preprocessing version of the xy bootstrap
#	globs as in "pfn" method
#	matrix s contains multinomial resampling draws (n by R)  
boot.rq.pxy <- function (x, y, s, tau = 0.5, coef, method="fn", Mm.factor = 3) {
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  R <- ncol(s)
  mofn <- nrow(s)
  #I have implemented the bootstrap with the count vector instead of the permutation vector
  counts <- apply(s, 2, tabulate, nbins=n)
  m <- round((p*n)^(1/2)) # Is there any theory for this?
  xxinv <- solve(chol(crossprod(x)))
  band <- sqrt(rowSums((x %*% xxinv)^2)) # H[i,i] vector
  r <- (y - x %*% coef)/band
  #To avoid computing quantiles of r in each replication I calculate once the cdf
  rs <- sort(r)
  B <- matrix(NA, R, p)
  #add two observations for the globs 
  x <- rbind(0,0,x)
  y <- c(0,0,y)
  counts <- rbind(1,1,counts)
  for(i in 1:R){
    mm <- m
    not.optimal <- TRUE
    #keep only the observations with positive weights
    bs <- (counts[,i] > 0)
    xb <- x[bs,]
    yb <- y[bs]
    rb <- r[bs[3:length(bs)]]
    wb <- counts[bs, i]
    while (not.optimal) {
      M <- Mm.factor * mm
      lo.q <- max(1, ceiling(n*tau - M/2))
      hi.q <- min(floor(n*tau + M/2), n)
      kappa <- rs[c(lo.q, hi.q)]
      sl <- c(FALSE, FALSE, rb < kappa[1])
      su <- c(FALSE, FALSE, rb > kappa[2])
      while (not.optimal) {
        if (any(sl)) {
          xb[1,] <- colSums(xb[sl, , drop = FALSE] * wb[sl])
          yb[1] <- sum(yb[sl] * wb[sl])
	} else sl[1] <- TRUE
        if (any(su)) {
          xb[2,] <- colSums(xb[su, , drop = FALSE] * wb[su])
          yb[2] <- sum(yb[su] * wb[su])
        } else su[2] <- TRUE
        sel <- c(!su & !sl)
        #use weighted QR
        z <- tryCatch(z <- rq.wfit(xb[sel, ], yb[sel], weights=wb[sel], 
		  tau = tau, method = method), error = function(e) e)
        if(inherits(z, "error")){
          mm <- 2*mm
          warning("Singular design: doubling m")
          break
          }        
        b <- z$coef
        rb <- yb[3:sum(bs)] - xb[3:sum(bs), ] %*% b
        su.bad <- c(FALSE, FALSE, rb < 0) & su
        sl.bad <- c(FALSE, FALSE, rb > 0) & sl
        bad.signs <- sum(su.bad | sl.bad)
        if (bad.signs > 0) {
          if (bad.signs > 0.1 * M) {
            mm <- 2*mm
            #warning("Too many fixups:  doubling m")
            break
          }
          su <- su & !su.bad
          sl <- sl & !sl.bad
        }
        else not.optimal <- FALSE
      }
    }
    B[i,] <- b
  }
  B
}

boot.rq.pwxy <- function (x, y, tau, coef, R = 200, m0 = NULL, eps = 1e-06, ...) 
{
    n <- length(y)
    if (nrow(x) != n) 
        stop("x and y don't match n")
    p <- ncol(x)
    if (!length(m0)) 
        m0 <- round(sqrt(n*p))
    r <- y - x %*% coef
    b <- matrix(0, p, R)
    nit <- matrix(0, 5, R) # This version keeps track of fixups!
    info <- rep(0, R)
    xxinv <- solve(chol(crossprod(x)))
    band <- pmax(eps, sqrt(((x %*% xxinv)^2) %*% rep(1, p)))
    r <- r/band
    o <- order(r)
    r <- r[o]
    y <- y[o]
    x <- x[o,]
    loq <- max(1, ceiling(n*tau - m0/2))
    hiq <- min(floor(n*tau + m0/2), n)
    qk <- r[c(loq,hiq)]
    z <- .Fortran("pwxy", as.integer(n), as.integer(p), as.integer(R), 
        as.double(t(x)), as.double(y), as.double(tau), as.double(qk),as.double(r), 
        b = as.double(-b), double(n), as.double(band), as.integer(m0), double(n), 
        double(n), double(n * 9), double(p * (p + 3)), double(p * 
            n), double(n), integer(n), integer(n), double(p), 
        double(p), double(p), nit = as.integer(nit), info = as.integer(info))
    if(any(z$info != 0)) warnings("Some bootstrap replications have abnormal convergence")
    B <- t(matrix(-z$b, p, R))
    nit <- matrix(z$nit,5,R)
    list(B = B, nit = nit)
}
