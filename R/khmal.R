"khmaladze.test" <-
function(fit, nullH = "location-scale", trim = c(0.25, 0.75) ) 
{
    # Compute Khmaladze test statistics:  fit is produced by rqProcess()
    J <- standardize(fit, nullH)
    Vtilde <- khmaladzize(fit$taus, fit$coef[1,], J$Vhat, nullH)
    vtilde <- khmaladzize(fit$taus, fit$coef[1,], J$vhat, nullH)
    trim <- ( fit$taus >= trim[1] & fit$taus <= trim[2] )	
    if(nullH == "location-scale") {
        Tvtilde <- (vtilde-vtilde[,2])/sqrt(max(fit$taus)-min(fit$taus))
        Tn <- max( apply(abs(Vtilde-Vtilde[,2])/
		 sqrt(max(fit$taus)-min(fit$taus)), 2,"sum")[trim] )
        THn <- apply(abs(Tvtilde[,trim]),1,max)
    	}
    else if( nullH ==  "location"){
        Tn <- max( apply( abs(Vtilde-Vtilde[,2])/
                sqrt(max(fit$taus)-min(fit$taus)), 2,"sum")[trim] )
        THn <- apply(abs((vtilde[,trim]-vtilde[,2])/
                sqrt(max(fit$taus)-min(fit$taus))),1,max)
        }
    x <- list(nullH = nullH, Tn = Tn, THn = THn)
    class(x) <- "khmaladze"
    x
}
"standardize" <-
function (fit, nullH = "location-scale") 
{
    Vhat <- fit$coefs
    vhat <- fit$coefs
    p <- nrow(vhat)
    if (nullH == "location-scale") {
	b <- matrix(0,2,p)
	ER <- fit$coefs
        for (i in 2:p) {
            er <- lm(fit$coefs[i, ] ~ fit$coefs[1, ])
            b[, i] <- er$coef
            ER[i, ] <- er$resid
	   }
        for (j in 1:length(fit$taus)) {
            V <- fit$Hinv[, , j] %*% fit$J %*% fit$Hinv[, , j]
            v <- V[-1, -1] + V[1, 1] * outer(b[2, -1], b[2, -1]) - 
                outer(V[-1, 1], b[2, -1]) - t(outer(V[-1, 1], b[2, -1]))
            v <- solve(v)
            v <- chol(0.5 * (v + t(v)))
            Vhat[-1, j] <- v %*% ER[-1, j]
            for (i in 2:dim(V)[1]) {
                v <- V[i, i] + V[1, 1] * b[2, i]^2 - 2 * V[i, 1] * b[2, i]
                vhat[i, j] <- ER[i, j]/sqrt(v)
                }
            }
        }
    else if (nullH == "location") {
        b <- apply(fit$coefs, 1, mean)
        ER <- fit$coefs - b
        for (j in 1:length(fit$taus)) {
            V <- fit$Hinv[, , j] %*% fit$J %*% fit$Hinv[, , j]
            v <- V[-1, -1]
            v <- solve(v)
            v <- chol(0.5 * (v + t(v)))
            Vhat[-1, j] <- v %*% ER[-1, j]
            for (i in 2:dim(fit$Hinv)[1]) {
                v <- V[i, i]
                vhat[i, j] <- ER[i, j]/sqrt(v)
                }
            }
        }
    else stop(paste("unrecognized nullH: ",nullH)) 
    list(Vhat=Vhat, vhat=vhat)
}
"khmaladzize" <-
function(tau, atau, Z, nullH)
{
	p <- diff( tau )
	p <- c( p[1], p )
	score <- akj(atau, atau, p)  
	L <- length(tau)
	gdot2 <- -score$psi
	gdot <- cbind(rep(1,L),gdot2)
	if(nullH == "location-scale")
	{
		gdot3 <- gdot2 * atau
		gdot <- cbind(gdot, gdot3)
	}
	kmin <- 0
	p <- nrow(Z)
	# the v process 
	for (i in 1:p) 
	{
		v <- Z[i,]
		dtau <- diff(tau)
		dtau <- c(dtau[1], dtau)
		dv <- c(0,diff(v))
		x1 <- gdot*sqrt(dtau)
		x1 <- x1[L:1,]
		y1 <- rev(dv/sqrt(dtau))
		bhat <- lm.fit.recursive(x1,y1,int=FALSE)
		bhat <- bhat[,L:1]
		dvhat <- diag(gdot%*%bhat)*dtau
		vhat <- cumsum(dvhat)
		v <- v[kmin:L-kmin]
		vhat <- vhat[kmin:L-kmin]
		Z[i,] <- v-vhat
	}
	return(Z)
}
"rqProcess" <-
function(formula, data, taus=seq(0.2,0.8,by=0.1) ) 
{
	z <- summary(lm(formula, data=data))
	ols <- z$coef
	vars <- names(z$coef[,1])
	p <- length(z$coefficients[,1])
	J <- solve(z$cov.unscaled)
	#
	# Compute RQ process 
	#
	ntaus <- length( taus )
	coefs<- matrix( 0, p, ntaus )
	Hinv <- array( 0, c(p, p, ntaus) )
	cat("taus: ")
	for(i in 1:ntaus)
	{
		cat(taus[i]," ")
		z <- summary(rq(formula,tau=taus[i],method="fn", data=data), se = "nid", hs=FALSE)
		coefs[,i] <- z$coefficients[,1]
		Hinv[,,i] <- z$Hinv
	}
    dimnames(coefs) <- list(vars,NULL)
    x <- list(taus = taus, coefs=coefs, J=J, Hinv=Hinv)
    class(x) <- "rqProcess"
    x
}
