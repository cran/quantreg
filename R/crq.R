crq.fit.pow <- function(x, y, yc, tau=0.5, weights=NULL, start, left=TRUE,maxit = 500){

 x <- as.matrix(x)
 y <- as.vector(y)
 n <- nrow(x)
 p <- ncol(x)
 if(missing(yc)) {
     yc <- rep(0,n) 
     if(!left) warning("default right (!?) censoring at zero")
    } 
 yc <- as.vector(yc)
 if(length(weights)){
     if (any(weights < 0))
         stop("negative weights not allowed")
     contr <- attr(x, "contrasts")
     x <- x * weights
     y <- y * weights
     }
 if(left) {
	x <- -x
	y <- -y
	yc <- -yc
	tau <- 1-tau
	}
# Starting value of beta 
if(missing(start)){ 
   tol <- .Machine$double.eps^(2/3)
   r <- rq.fit.br(x,y,tau)$residuals
   h <- which(abs(r) < tol)[1:p]
   }
else if(start[1] == "global"){
# Generate H matrix 
   H <- combos(n,p)
   m <- ncol(H)

# Starting value of beta 
   U <- solve(x[H[,1],])
   bh <- U %*% y[H[,1]]

f <- .Fortran("brutpow",
        as.integer(n),
        as.integer(p),
        as.integer(m),
        as.integer(H),
        as.double(x),
        as.double(y),
        as.double(yc),
        as.double(bh),
        as.double(tau),
        as.double(U),
        double(p),
        double(p),
        kminz = integer(1),
        nflag = as.integer(0))
if(f$nflag!=0)
        warning(switch(f$nflag,
                "Error in pivot:  hout not in h",
                "Error in pivot:  hin in h",
                "Error in pivot:  hin out of bounds",
                "Error in findk:  k not found"
                ))
   k  <- f$kminz
   bh <- solve(x[H[,k],],y[H[,k]])
   residuals <- as.matrix(y - pmin(yc, x %*% bh))

   return(list(coefficients= bh,residuals = residuals))
  }
else if(is.numeric(start) && length(start) == p){
	if(all(start %in% 1:n)) h <- start
	else  h <- order(abs(y - x %*% start))[1:p]
	}
else
   stop("Invalid starting value")

xhinv <- try(solve(x[h,]))
if(class(xhinv)[1] == "try-error")
	stop("Singular basic solution generated by 'start'")

f <- .Fortran("powell", 
	as.integer(n), 
	as.integer(p), 
	as.integer(2*p), 
	as.double(x),
	as.double(y), 
	as.double(yc), 
	coef = double(p), 
	as.double(tau), 
	as.integer(h), 
	double(n), 
	as.double(xhinv), 
	double(n), 
	double(2*p),
	double(p),
	double(p),
	as.integer(maxit), 
	nflag = as.integer(0))
	
if(f$nflag!=0) 
	warning(switch(f$nflag, "Max iterations reached",
		"Solution may be nonunique",
		"Error in pivot:  hout not in h",
		"Error in pivot:  hin in h",
		"Error in pivot:  hin out of bounds"
		))
coef <- f$coef
residuals <- as.matrix(y - pmin(yc, x %*% coef))

return (list(coefficients= coef, residuals = residuals))
}


Curv <- function (y,  yc, ctype = c("left", "right")) 
{
    nn <- length(y)
    if(length(yc) != nn)
        stop("Event times and censoring times of different length")
    ctype <- match.arg(ctype)
    if(ctype == "right" && any(y > yc))
        stop("Event times can not exceed ctimes for right censoring")
    if(ctype == "left" && any(y < yc))
        stop("Event times can not be less than ctimes for left censoring")
    if(!(ctype == "left" || ctype == "right"))
        stop("Invalid ctype for method Powell")
    ss <- cbind(y, yc)
    dimnames(ss) <- list(NULL, c("time", "ctime"))
    attr(ss, "ctype") <- ctype
    attr(ss, "type") <- ctype
    class(ss) <- "Surv"
    ss
}

QTECox <- function(x, smooth = TRUE){
# compute quantile treatment effect for a Cox PengHuang fit
g <- survival::survfit(x)
if(smooth)
        g <- supsmu(1-g$surv,g$time)
taus <- (g$x[-1] + g$x[-length(g$x)])/2
Qhat <- (g$y[-1] + g$y[-length(g$y)])/2

dQ <- diff(Qhat)/diff(taus)
taus <- taus[-1]; Qhat <- Qhat[-1]
QTE <- outer(dQ * (1-taus) * log(1-taus) / Qhat, coef(x))
list(QTE = QTE , taus = taus)
}
boot.crq <- function(x, y, c, taus, method, ctype = "right", R=100,  
		mboot,  bmethod = "jack", ...)
        {
        n <- length(y)
        p <- ncol(x)
        if(missing(mboot)) {
                if(bmethod=="jack") mboot <- 2*ceiling(sqrt(n))
                else mboot <- n }
	if(length(taus) > 1)
	    A <- array(0,dim=c(p,length(taus),R))
	else
	    A <- matrix(0,p,R)
        for (i in 1:R){
                if(bmethod == "jack") { 
                        s <- sample(1:n,mboot)
                        yb <- y[-s]
                        xb <- x[-s, ]
                        cb <- c[-s]  
                        w <- rep(1,n-mboot)
                        }
                else if(bmethod == "xy-pair"){
                        w <- table(sample(1:n,mboot,replace=TRUE))
                        s <- as.numeric(names(w))
                        w <- as.numeric(w)
                        yb <- y[s]
                        xb <- x[s,]
                        cb <- c[s]
                        }
                else if(bmethod == "Bose"){
                        w <- rexp(n)
                        yb <- y
                        xb <- x
                        cb <- c
                        }
                else
                        stop("invalid bmethod for boot.crq")
		if(method == "Portnoy")
                	a <- crq.fit.por(xb,yb,cb, weights = w, ctype = ctype, ... )
		#else if(method == "Portnoy2")
                	#a <- crq.fit.por2(xb,yb,cb, weights = w, ctype = ctype, ... )
		else if(method == "PengHuang")
                	a <- crq.fit.pen(xb,yb,cb, weights = w, ctype = ctype, ... )
		else
			stop("Invalid method for boot.crq")
                if((i %% floor(R/10)) == 0 & n > 100000)
                        cat(paste("bootstrap roughly ",100*(i/R)," percent complete\n"))
		if(length(taus) > 1)
		    A[,,i] <- coef(a,taus)
		else
		    A[,i] <- coef(a,taus)
                }
        list(A = A, n = length(y), mboot = mboot, bmethod = bmethod)
        }

coef.crq <- function(object, taus = 1:4/5, ...)
        {
        # Extract coefficients from the crq solution array 
        
        if(min(taus) < 0 || max(taus) > 1) stop("taus out of range [0,1]")
	if(length(object$tau) == 1){
        	coef <- object$coefficients
        	return(coef)
        	}
        taus <- sort(taus)
        S <- object$sol
        ctype <- object$ctype
        r <- S[1, ]
        r <- c(r[1],r)
	if(is.unsorted(r)) r <- r[-length(r)] #kludge until why this happens is found
        B <- S[-1,,drop = FALSE]
	B <- t(cbind(B,B[,ncol(B),drop = FALSE]))
        ts <- taus[taus > min(r) & taus < max(r)]
        bin <- findInterval(ts,r)
        wgt <- (ts - r[bin])/(r[bin + 1] - r[bin])
        binlag <- bin - 1
        binlag[binlag == 0] <- 1
        coef <- t(wgt * B[bin, , drop = FALSE] + (1 - wgt) * B[binlag, , drop = FALSE])
        nna <- length(taus) - length(ts)
        if(nna > 0) {
           if (ctype == "left")
              coef <- cbind(matrix(NA, nrow(coef), nna), coef)
           else
              coef <- cbind(coef, matrix(NA, nrow(coef), nna))
        }
        taulabs <- paste("tau=", format(round(taus, 3)))
        dimnames(coef)[[2]] <- taulabs
        coef[-nrow(coef),]  # Delete Qhat entries
        }

crq <- function (formula, taus, data, subset, weights, na.action, 
	method = c("Powell","Portnoy","Portnoy2", "PengHuang"), contrasts = NULL, ...)
{
    if(!requireNamespace("survival", quietly = TRUE))
	warning("crq requires survival package to be installed")
    #if(method == "Portnoy2") stop("Portnoy2 method not (yet) implemented")
    Surv <- survival::Surv
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
        names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    if (method == "model.frame")
        return(mf)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf, contrasts)
    weights <- as.vector(model.weights(mf))
    Y <- model.extract(mf, "response")
    eps <- .Machine$double.eps^(2/3)
    if(!inherits(Y,"Surv"))
        stop("Response must be a survival object")
    method <- match.arg(method)
    if(method == "Powell") {
	ctype <- attr(Y, "ctype")
	if(!ctype %in% c("right","left")) 
             stop("Only right or left censoring Surv objects are allowed")
	left <- (ctype == "left")
        if (any(taus < -eps) || any(taus > 1 + eps))
            stop("invalid taus:  taus should be >= 0 and <= 1")
	y <- Y[,1]
	cen <- Y[,2]
	if (length(taus) > 1) {
          coef <- matrix(0, ncol(X), length(taus))
          fitted <- resid <- matrix(0, nrow(X), length(taus))
          for (i in 1:length(taus)) {
            z <- crq.fit.pow(X, y, cen, tau = taus[i], weights, left = left, ...)
            coef[, i] <- z$coefficients
            resid[, i] <- z$residuals
            fitted[, i] <- y - z$residuals
          }
          taulabs <- paste("tau=", format(round(taus, 3)))
          dimnames(coef) <- list(dimnames(X)[[2]], taulabs)
          dimnames(resid) <- list(dimnames(X)[[1]], taulabs)
          fit <- list(coefficients = coef, residuals = resid, fitted.values = fitted)
	  fit$tau <- taus
	  class(fit) <- "crqs"
	}
	else {
          fit <- crq.fit.pow(X, y, cen, tau = taus, weights, left = left, ...)
	  fit$tau <- taus
	  class(fit) <- "crq"
	}
     }
    else if(method == "Portnoy"){
        ctype <- "right"
        if(attr(Y,"type") != "right"){
            if(attr(Y,"type") == "left") ctype <- "left"
            else stop("Only right censoring Surv objects are allowed for Portnoy method")
            }
        y <- Y[,1]
        cen <-  Y[,2]
        fit <- crq.fit.por(X, y, cen, weights, ctype = ctype,   ...)
	class(fit) <- "crq"
	}
    else if(method == "Portnoy2"){
        ctype <- "right"
        if(attr(Y,"type") != "right"){
            if(attr(Y,"type") == "left") ctype <- "left"
            else stop("Only right censoring Surv objects are allowed for Portnoy method")
            }
        y <- Y[,1]
        cen <-  Y[,2]
        fit <- crq.fit.por2(X, y, cen, weights, ctype = ctype,   ...)
	class(fit) <- "crq"
	}
    else if(method == "PengHuang"){
        ctype <- "right"
        if(attr(Y,"type") != "right"){
            if(attr(Y,"type") == "left") ctype <- "left"
            else stop("Only right censoring Surv objects are allowed for Portnoy method")
            }
        y <- Y[,1]
        cen <-  Y[,2]
        fit <- crq.fit.pen(X, y, cen, weights,  ctype = ctype, ...)
	class(fit) <- "crq"
	}
    else
	stop("Method not defined for crq")

fit$terms <- mt
fit$call <- call
fit$formula <-  formula(mt)
fit$method <-  method
fit$contrasts <- contrasts 
fit$ctype <-  ctype
attr(fit, "na.message") <- attr(m, "na.message")
fit
}
predict.crqs <-
function (object, newdata, type = NULL, ...) {
pred <-  predict.rqs(object, newdata, ...) 
	}
predict.crq <-
function (object, newdata,  ...) {
method <- object$method
if(method %in% c("Portnoy","Portnoy", "PengHuang")){ #Kludge to make crq sol matrix look like rq sol matrix 
	p2 <- nrow(object$sol)
	object$sol <- object$sol[c(1,p2,p2,2:(p2-1)),]
	}
pred <- switch(method,
	"Powell" = predict.rq(object, newdata,  ...), 
	"Portnoy"  = predict.rq.process(object, newdata,  ...), 
	"Portnoy2"  = predict.rq.process(object, newdata,  ...), 
	"PengHuang" = predict.rq.process(object, newdata,  ...) 
	)
}
crq.fit.por <- function(x, y, cen, weights = NULL, grid, ctype = "right")
{
      if(!is.matrix(x))
	  stop("x must be a matrix")
      if(!length(dimnames(x)))
	  dimnames(x)[[2]] <- paste("V",1:p)
      p <- ncol(x)
      n <- length(y) 
      cen <- 1 - cen #NB: Fortran routine wants censoring indicator flipped (!!!)
      mp <- max(n + 5 + max(1, sum(cen)),n+p)
      eps <- 1e-04
      kmax <- 10000 # perhaps this should be parameter in future
      if(length(weights)){
		if (any(weights < 0)) 
		        stop("negative weights not allowed")
    		contr <- attr(x, "contrasts")
    		x <- weights * x
    		y <- weights * y
		}
    if(ctype == "left") y <- -y
    if(missing(grid))  grid <-  seq(1/n,1-1/n,by=1/(5 + 3 * n^.4))
    if(is.numeric(grid)){
		ginit <- min(grid)
		dgrid <- diff(grid)
		gstep <- median(dgrid)
		if(any(dgrid <0)) stop("grid is not monotonic")
		if(gstep < eps) stop("grid stepsize too small")
		nsol <- 3*n
		mw = -1
		}
      else if(grid  == "pivot"){
		nsol <- 3*n
		ginit <- 1/(2*n)
		gstep <- 1/(2*n)
		mw <- min(10, mp) # Somewhat arbitrary see crq.f
		}
      else 
		stop("Invalid grid")
	z <- .Fortran("crqf",
		as.integer(n),
		as.integer(p),
                as.integer(mp), 
		as.integer(p+2),
		as.double(x),
		as.double(y),
		as.integer(cen),
		as.double(ginit),
                as.integer(mw),
                as.double(gstep),
		ift = integer(1),
		h = integer(p),
                xh = double(p*p),
		wa = double(mp*p),
		wb = double(mp),
		wc = double(mp*(p+2)),
		wd = double(mp),
		we = double(mp),
		wf = double(p),
		iflag = integer(mp),
                as.integer(nsol),
		sol = double(nsol*(p+2)),
                lsol = integer(1),
		icen = integer(n),
		tcen = double(n),
		lcen = integer(1))
	nw <- z$h[1]
	flag <- z$ift
	msg <- switch(flag,
		paste("Error in input dimensions, n,p,mw "),
		paste("Error in input dimensions, n,p,mw "),
		paste("Error in input dimensions, n,p,mw "),
		paste("Less than p=",p,"observations above tau = 0 solution"),
		paste("Possible degeneracy at",nw,"tau values.",
          		"$tau.degen: first mp =", n + 5 + sum(cen)," such tau values"), 
		paste("Number of pivots to be saved in sol > nsol.",
          		"Redefine nsol: use nsol < n to save for tau = i/(nsol-1)"),
		paste("Error with partial return: possible degeneracies",
         		"Max number of rq calls exceeded: dither x or increase mw"),
		paste("Premature stop: defective conditional distribution"),
		paste("Simplex iteration limit exceeded -- consider dithering y"))
	#if(flag > 0 && flag != 5 && flag < 10)
	if(flag %in% c(1:4,6,7))
		ifelse(flag <= 3,stop(msg),warning(msg))
	J <- z$lsol
	B <- matrix(z$sol, nrow=p+2, ncol=nsol, byrow=FALSE)[,1:J, drop = FALSE]
	if(B[1,J] < B[1,J-1]) B <- B[,-J] # SLP hack Oct 27 2014
	ic <- z$icen
	sp <- (1:n)[ic == 1]
	tsp <- z$tcen[sp]
        t1 <- z$wd[1:nw]
        if(ctype == "left") {
            B[1,] <- 1 - B[1,]
            B[-1,] <- - B[-1,]
            B <- B[,ncol(B):1]
            }
	dimnames(B) <- list(c("tau",dimnames(x)[[2]],"Qbar"),NULL)
	a <- list(sol=B, Isplit=sp, tsp = tsp, status=ic, ctype = ctype)
	fitted <- x %*% B[-c(1,dim(B)[1]),]
	dimnames(fitted) <- list(NULL, paste("tau=",round(B[1,],4))) 
	a$fitted.values <- fitted
	class(a) <- "crq"
	return(a) 
	}
crq.fit.pen <- function(x, y, cen, weights=NULL,grid, ctype = "right" ){
      p <- ncol(x)
      n <- length(y)
      if(missing(grid)) grid <-  seq(1/n,1-1/n,by=min(0.01,1/(2*length(y)^.7))) 
      if(!is.numeric(grid)) stop("Invalid grid")
      if(any(grid < 0) || any(grid > 1)) stop("Invalid grid")
      m <- length(grid)
      xbar <- apply(x,2,mean)
      if(length(weights)){
                if (any(weights < 0))
                        stop("negative weights not allowed")
                contr <- attr(x, "contrasts")
                x <- x * weights
                y <- y * weights
                }
        if(ctype == "left") y <- -y
        s <- rep(0,n)
        u <- rep(1,n)
        d <- rep(1,n)
        r <- rep(1,p)
        B <- matrix(0,p,m)
        cc <- as.logical(cen)
        y1 <- y[cc]
	n1 <- length(y1)
        x1 <- x[cc,]
	z <- .Fortran("crqfnb", as.integer(n), as.integer(p), 
		a1 = as.double(t(as.matrix(x1))), c1 = as.double(-y1), n1=as.integer(n1),
		as.double(x), as.double(y),as.double(cen),B =as.double(B),
		g = as.double(grid),m = as.integer(m), as.double(r), 
		as.double(s), as.double(d), as.double(u),
        	wn = double(n1 * 9), wp = double((p + 3) * p),
        	info = integer(1))
	J <- z$m - 1
	B <- matrix(-z$B, p, m)
	B <- B[,1:J,drop = FALSE]
        qhat <- t(xbar) %*% B
        B <- rbind(grid[1:J],B,qhat)
        dimnames(B) <- list(c("tau",dimnames(x)[[2]],"Qhat"),NULL)
        if(ctype == "left") {
            B[1,] <- 1 - B[1,]
            B[-1,] <- - B[-1,]
            B <- B[,ncol(B):1]
            }
        B  <- list(sol=B, ctype = ctype)
	class(B) <- "crq"
	B
        }
plot.summary.crqs <-
function (x, nrow = 3, ncol = 3, CoxPHit = NULL, ...) {
    taus <- function(x) x$tau
    xx <- unlist(lapply(x, taus))
    coef <- lapply(x, coefficients)
    dots <- list(...)
    p <- nrow(coef[[1]])
    k <- ncol(coef[[1]])
    if(k != 6) stop("summary.crqs object has wrong column dimension")
    m <- length(xx)
    blab <- dimnames(coef[[1]])[[1]]
    a <- array(unlist(coef), c(p, k, m))
    if(length(CoxPHit))
	CoxQTE <- QTECox(CoxPHit)
    oldpar <- par(no.readonly=TRUE)
    par(mfrow = c(nrow, ncol))
    for (i in 2:p) {
         b  <- a[i, 1, ]
         bl <- a[i, 2, ]
         bu <- a[i, 3, ]
	if(length(dots))
	    plot(rep(xx, 2), c(bl, bu), type = "n", ...)
	else
	    plot(rep(xx, 2), c(bl, bu), xlab = "", ylab = "", type = "n")
        title(paste(blab[i]), cex = 0.75)
        polygon(c(xx, rev(xx)), c(bl, rev(bu)), col = "LightSkyBlue")
        points(xx, b, cex = 0.5, pch = "o", col = "blue")
        lines(xx, b, col = "blue")
        abline(h = 0)
        if(length(CoxPHit)) {
	    lines(CoxQTE$taus,CoxQTE$QTE[,i-1],col="red")
            }
    }
par(oldpar)
}
summary.crq <-
function (object, taus = 1:4/5, alpha = .05, se = "boot", covariance = TRUE, ...) {
    mt <- terms(object)
    m <- model.frame(object)
    x <- model.matrix(mt, m, contrasts = object$contrasts)
    Y <- model.response(m)
    method <- object$method
    ctype <- object$ctype
    y <- Y[,1]
    cen  <- Y[,2]
    wt <- as.vector(model.weights(object$model))
    if (!is.null(wt)) {
        resid <- resid * wt
        x <- x * wt
        y <- y * wt
	cen <- cen * wt
    }
    if(method == "Powell"){
       coef <- coefficients(object)
       vnames <- dimnames(x)[[2]]
       resid <- object$residuals
       tau <- object$tau
       n <- length(resid)
       p <- length(coef)
       rdf <- n - p
       if(rdf < 1) stop("residual degrees of freedom nonpositive")
       if (se == "boot") {
           if(attr(Y,"type") == "left")
              s <- cen  <  x %*% coef
           else
              s <- cen  >  x %*% coef
           B <- boot.rq(x[s, ], y[s], tau, ...)
           cov <- cov(B$B)
           serr <- sqrt(diag(cov))
           }
       else stop("Only boot method is implemented for crq inference")
       coef <- array(coef, c(p, 6))
       fact <- qt(1 - alpha/2, rdf)
       coef[, 2] <- coef[,2] - fact * serr
       coef[, 3] <- coef[,3] + fact * serr
       coef[, 4] <- serr
       coef[, 5] <- coef[, 1]/coef[, 4]
       coef[, 6] <- 2 * (1 - pt(abs(coef[, 5]),rdf)) 
       cnames <- c("Value","Lower Bd","Upper Bd","Std Error","T Value","Pr(>|t|)")
       dimnames(coef) <- list(vnames, cnames)
       object <- object[c("call", "terms")]
       if (covariance == TRUE) 
           object$cov <- cov
       object$B <- B
       object$coefficients <- coef
       object$rdf <- rdf
       object$tau <- tau
       class(object) <- "summary.crq"
       return(object)
       }
    else if(method == "Portnoy" || method == "Portnoy2" || method == "PengHuang") {
	if(length(taus) == 1) taus <- c(taus, taus) # Kludge to avoid dimensional errors
       coef <- as.matrix(coef(object,taus))
       coef <- coef[,apply(coef,2,function(x) any(!is.na(x))),drop = FALSE] # Delete NA columns if any
       if(ctype == "right") taus <- taus[1:ncol(coef)]
       else  taus <- taus[(1 + length(taus)-ncol(coef)):length(taus)]
       B <- boot.crq(x, y, cen, taus, method = method, ctype = ctype, ...)
       bmethod <- B$bmethod
       nas <- apply(is.na(B$A[1,,, drop = TRUE]),1,sum)
       Brep <- dim(B$A)[3]
       p <- dim(B$A)[1]
       if(bmethod == "jack") 
	   sqmn <- sqrt((B$n - B$mboot - p + 1)/B$mboot)
       else 
	   sqmn <- sqrt(B$mboot/B$n)
       fact <-   qnorm(1 - alpha/2)/qnorm(.75)
       A <- apply(B$A, 1:2, quantile, probs = 1:3/4, na.rm = TRUE)
       #D <- .5 * fact *(A[3,,]-A[1,,]) * sqmn
       #L <- coef - D
       #U <- coef + D
       DU <- fact *(A[3,,]-A[2,,]) * sqmn
       DL <- fact *(A[2,,]-A[1,,]) * sqmn
       L <- coef - DL
       U <- coef + DU
       S <- (U - L)/(2 * qnorm(1 - alpha/2))
       T <- coef/S
       P <- 2 * (1 - pnorm(abs(T)))
       G <- list()
       cnames <- c("Value","Lower Bd","Upper Bd","Std Error","T Value","Pr(>|t|)")
       for(i in 1:length(taus)){
	   tab <- cbind(coef[,i],L[,i],U[,i],S[,i],T[,i],P[,i])
	   dimnames(tab)[[2]] <- cnames
	   if(covariance) 
	       cov <- var(t(B$A[,i,]), na.rm = TRUE) * sqmn^2
	   else
	       cov <- NULL
	   G[[i]] <- list(tau = taus[i], coefficients = tab, NAs = nas[i], 
			  cov = cov, Brep = Brep, bmethod = bmethod)
	   }
       class(G) <- "summary.crqs"
       return(G)
       }
    else
       stop("Invalid method for summary.crq")
   }
"summary.crqs" <-
function (object, ...) {
	if(!object$method == "Powell")
		stop("Invalid method")
        taus <- object$tau
        xsum <- as.list(taus)
        for(i in 1:length(taus)){
                xi <- object
                xi$coefficients <- xi$coefficients[,i]
                xi$residuals <- xi$residuals[,i]
                xi$tau <- xi$tau[i]
                class(xi) <- "crq"
                xsum[[i]] <- summary(xi, ...)
                }
        class(xsum) <- "summary.crqs"
        xsum
        }

print.crq <- function(x, ...){
	if(!is.null(cl <- x$call)) {
                cat("Call:\n")
                dput(cl)
		}
        coef <- coef(x, ...)
        cat("\nCoefficients:\n")
	print.default(coef)
        invisible(x)
	}

print.summary.crqs <- function(x, ...)
    lapply(x,print.summary.crq)
print.summary.crq <- function (x, digits = max(5, .Options$digits - 2), ...) {
    coef <- x$coefficients
    tau <- x$tau
    if(length(x$NAs)) NAs <- x$NAs
    else NAs <- 0
    cat("\ntau: ")
    print(format(round(tau, digits = digits)), quote = FALSE, ...)
    if(NAs > 0) 
       cat(paste("   Number of NA Bootstrap Replications: ", NAs, "out of", x$Brep))
    cat("\nCoefficients:\n")
    print(format(round(coef, digits = digits)), quote = FALSE, ...)
    invisible(x)
}

#  SLP: Mar 31 2015 (Revised RWK:  Aug  2015)
#  An R version of crq.fit.por without any frills on a fixed grid.
crq.fit.por2 <- function(x, y, cen, weights = NULL, grid, ctype = "right") {
    method <- "fn"
    BIG <- 1e+6
    ztol <- 1e-6
    n <- length(y)
    p <- ncol(x)
    cen <- as.logical(cen)
    uncen <- !cen
    if(sum(uncen) < p)  
	stop("too few uncensored observations") 
    sol <- matrix(NA, p+2,1)
    if (length(weights)) {
        if (any(weights < 0)) 
            stop("negative weights not allowed")
        contr <- attr(x, "contrasts")
        x <- weights * x
        y <- weights * y
    }
    w <- rep(1,n)
    if (ctype == "left")  y <- -y
    if(missing(grid)) { # Fixme:  else if ...
	g0 <- 1/(2*n)
	dg <- 1/(5 + 3*n^.4)
	grid <- seq(g0,1-g0,by=dg)
    }
    repeat { # Peel censored points below g0 fit
	z <- rq.wfit(x, y, tau = g0, weights = w, method = method)
	J <- which(abs(z$resid) < ztol & cen)
	if(length(J)){
	    y <- y[-J]
	    x <- x[-J,]
	    w <- w[-J]
	    cen <- cen[-J]
	}
	else break
    }
    n <- length(y)
    U <- which(cen)
    if(length(cen) == 0) { # No censored points left!
	z <- rq(y ~ x - 1, tau = grid, weights = w, method = method)
	return(coef(z))
    }
    else { # Data Augmentation:
	m <- sum(cen)
	x <- rbind(x, x[cen, ])
	y <- c(y, rep(BIG, m))
	w[cen] <- 0
	w <- c(w, rep(1, m))
    }
    taw <- rep(0, length(y))
    K <- NULL
    sol <- matrix(NA, p+2, length(grid))
    b <- z$coef
    xbar <- apply(x, 2, mean)
    sol[,1] <- c(0, b, sum(xbar * b))
    for(i in 2:length(grid)){
	tau <- grid[i]
	if(length(K)){ # Update w
	    w[K] <- (tau - taw[K])/(1 - taw[K])
	    idx <- which(which(cen) %in% K)[rank(K)]
	    w[n + idx] <- 1 - w[K]
	}
	z <- rq.wfit(x, y, tau = tau, weights = w, method = method)
	b0 <- b
	b <- z$coef
	sol[,i] <- c(tau, b, sum(b * xbar))
	if(any(abs(z$resid) < ztol & (y > BIG - 1))) break
	J <- U[z$resid[U] < 0] 
	if(length(J)){
	    U <- setdiff(U, J)
	    K <- union(K, J)
	    y0 <- x[J,,drop = FALSE] %*% b0
	    y1 <- y[J] - z$resid[J]
	    gstep <- grid[i] - grid[i-1]
	    taw[J] <- tau - gstep * (y1 - y[J])/(y1 - y0)
	}
    }
    if (ctype == "left") {
        sol[1, ] <- 1 - sol[1, ]
        sol[-1, ] <- -sol[-1, ]
    }
    dimnames(sol) <- list(c("tau", dimnames(x)[[2]], "Qbar"), NULL)
    return(list(sol = sol, ctype = ctype)) 
}

