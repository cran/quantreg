"rqProcess" <-
function (formula, data, taus, nullH = "location", ...) 
{
    z <- summary(f <- lm(formula, data = data, x = TRUE))
    xbar <- apply(f$x,2,mean)
    vars <- names(z$coef[-1, 1])
    p <- length(z$coef[, 1])
    n <- length(z$resid)
    Jinv <- z$cov.unscaled
    pivot <- any(taus < 0) || any(taus > 1)
    if(!pivot){ #grid method
	if(abs(diff(range(diff(taus)))) > sqrt(.Machine$double.eps))
		stop("rqProcess must be evaluated on equally spaced taus")
	ntaus <- length(taus)
	z <- summary(rq(formula, data = data, tau = taus, method = "fn"), covariance = TRUE, ...)
	coefs = t(sapply(z, coefficients)[1:p,])
	Cov = array(sapply(z, function(x) x$cov), c(p,p,ntaus))
	qtaus <- coefs %*% xbar
	Vhat <- t(coefs)[-1,,drop=FALSE]
	vhat <- t(coefs)[-1,,drop=FALSE]
	J <- solve(Jinv)
	p <- nrow(J)
	if (nullH == "location-scale") {
	    f <- lm(coefs[,-1] ~ coefs[,1]) 
	    b <- matrix(f$coef,2,p-1)[2,]
	    R <- matrix(f$resid,ntaus,p-1)
	    for (j in 1:length(taus)) {
		V <- Cov[, , j]
		v <- V[-1, -1] + V[1, 1] * outer(b,b) - 
		    outer(V[-1, 1], b) - t(outer(V[-1, 1], b))
		v <- solve(v)
		v <- chol(0.5 * (v + t(v)))
		Vhat[,j] <- v %*% R[j,]
		for (i in 2:p) {
		    v <- V[i, i] + V[1, 1] * b[i-1]^2 - 2 * V[i, 1] * b[i-1]
		    vhat[i-1,j] <- R[j, i-1]/sqrt(v)
		    }
		}
	    }
	else if (nullH == "location") {
	    b <- apply(coefs, 2, mean)
	    R <- t(coefs) - b
	    for (j in 1:length(taus)) {
		V <- Cov[, , j] 
		A <- solve(V[-1, -1,drop=FALSE])
		B <- chol(0.5 * (A + t(A)))
		Vhat[,j] <- B %*% R[-1, j,drop=FALSE]
		vhat[,j] <- R[-1, j,drop=FALSE] / (sqrt(diag(V))[-1])
		}
	    }
	}
   else{
       stop("pivot method is now deprecated for Ktest")
	}
    dimnames(Vhat) <- list(vars, NULL)
    dimnames(vhat) <- list(vars, NULL)
    x <- list(taus = taus, qtaus = qtaus, Vhat = Vhat, vhat = vhat, n = n)
    class(x) <- "rqProcess"
    x
}
"KhmaladzeTest"  <-
function (formula, data = NULL, taus = 1:99/100, nullH = "location", 
	trim = c(0.05, 0.95), h = 1, ...) 
{
    gdot = function(x, taus = 1:999/1000, h = -1){
	Qt = quantile(x, taus)
	z = akj(x, Qt, h = h)
	Qf = approxfun(taus, Qt)
	den = approxfun(Qt, z$den)
	psi = approxfun(Qt, z$psi)
	gdot0 = den(Qf(taus))
	gdot1 = psi(Qf(taus))
	gdot2 = Qf(taus) * gdot1 + 1
	cbind(gdot0, gdot1, gdot2)
    }
    f <- rqProcess(formula, data = data, taus=taus, nullH = nullH, ...)
    G = gdot(f$qtaus, taus, h = 1)
    ntaus = length(taus)
    Vtil = Vhat = G[,1] * f$Vhat
    vtil = vhat = G[,1] * f$vhat
    if(nullH == "location") X = G[-ntaus,2] 
    else X = G[-ntaus,2:3]
    for(i in 1:nrow(Vhat)){
	Y = Vhat[i,-1]
	y = vhat[i,-1]
	B = lm.fit.recursive(X,Y)
	b = lm.fit.recursive(X,y)
	Vtil[i,-1] = Y - diag(cbind(1,X) %*% B)
	vtil[i,-1] = y - diag(cbind(1,X) %*% b)
    }
    trim <- (f$taus >= trim[1] & f$taus <= trim[2])
    Tvtil <- (vtil - vtil[, 2])/sqrt(max(f$taus) - min(f$taus))
    TVtil <- apply(abs(Vtil - Vtil[, 2])/
	sqrt(max(f$taus) - min(f$taus)), 2, "sum")[trim]
    Tn <-  max(TVtil)
    THn <- apply(abs(Tvtil[, trim,drop = FALSE]), 1, max)
    x <- list(nullH = nullH, Tn = Tn, THn = THn, taus = taus, 
	      Vhat = t(Vhat), Vtil = t(Vtil))
    class(x) <- "KhmaladzeTest"
    x
}
print.KhmaladzeTest = function(x, ...) {
    cat("\nTest of H_0:", x$nullH, "\n")
    cat("\nJoint Test Statistic:", x$Tn, "\n")
    cat("\nComponent Test Statistics:", x$THn, "\n")
}
plot.KhmaladzeTest = function(x, ...) {
    dev.new(height = 6, width = 10)
    par(mfrow = c(1,2))
    matplot(x$taus, x$Vhat, xlab = expression(tau), ylab = expression(hat(v)[n](tau)), type = "l")
    title("Parametric QR Process")
    matplot(x$taus, x$Vtil, xlab = expression(tau), ylab = expression(tilde(v)[n](tau)), type = "l")
    title("Transformed Parametric QR Process")
    par(mfrow = c(1,1))
}
