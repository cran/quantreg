"anova.rq" <-
function (object, ...)
{
    if (length(list(object, ...)) > 1) {
        return(anova.rqlist(object, ...))
    }
    stop("Anova is only defined (yet) for sequences of rq objects")
}
"anova.rqlist" <-
function (object, ..., test = "Wald") 
{
    objects <- list(object, ...)
    responses <- as.character(lapply(objects, function(x) formula(x)[[2]]))
    sameresp <- responses == responses[1]
    if (!all(sameresp)) 
        stop("Models don't all have the same response variable")
    n <- length(objects[[1]]$y)
    models <- as.character(lapply(objects, function(x) formula(x)))
    if (test == "Wald") 
        objects <- lapply(objects, function(x) summary(x))
    nobjects <- length(objects)
    dimp <- lapply(objects, function(x) length(x$coef[, 1]))
    objects <- objects[order(-unlist(dimp))]
    models <- as.character(lapply(objects, function(x) formula(x)))
    taus <- unlist(lapply(objects, function(x) x$tau))
    names <- lapply(objects, function(x) names(x$coefficients[, 
        1]))
    sametaus <- taus == taus[[1]]
    if (all(sametaus)) {
        Tn <- rep(0, nobjects - 1)
        ndf <- Tn
        ddf <- Tn
        pvalue <- Tn
        topnote <- paste("Model ", format(1:nobjects), ": ", 
            models, sep = "", collapse = "\n")
        if (test == "Hajek") {
            x1 <- as.matrix(objects[[1]]$x)
            y <- objects[[1]]$y
            for (i in 2:nobjects) {
                if (!all(names[[i]] %in% names[[1]])) 
                  stop("Models aren't nested")
                nullH <- is.na(match(names[[1]], names[[i]]))
                X1 <- as.matrix(x1[, nullH])
                X0 <- as.matrix(objects[[i]]$x)
                Htest <- rq.test.Hajek(X0, X1, y, score = "wilcoxon")
                Tn[i - 1] <- Htest$Tn
                ndf[i - 1] <- Htest$ndf
                ddf[i - 1] <- Htest$ddf
                pvalue[i - 1] <- Htest$pvalue
            }
        }
        else if (test == "Wald") {
            V <- lapply(objects, function(x) x$cov)
            coef <- lapply(objects, function(x) x$coefficients[, 
                1])
            for (i in 2:nobjects) {
                if (!all(names[[i]] %in% names[[1]])) 
                  stop("Models aren't nested")
                nullH <- is.na(match(names[[1]], names[[i]]))
                ndf[i - 1] <- sum(nullH)
                Tn[i - 1] <- t((coef[[1]])[nullH]) %*% solve((V[[1]])[nullH, 
                  nullH], (coef[[1]])[nullH])/ndf[i - 1]
                ddf[i - 1] <- n - length(names[[1]])
                pvalue[i - 1] <- 1 - pf(Tn[i - 1], ndf[i - 1], 
                  ddf[i - 1])
            }
        }
        else stop("Mode test only defined for Wald and Hajek")
    }
    else {
        m <- length(taus)
        for (i in 2:m) {
            if (!setequal(names[[i]], names[[1]])) 
                stop("Models with common tau don't have same X")
        }
        if (names[[1]] != "(Intercept)") 
            stop("Intercept required in common tau testing")
        Omega <- outer(taus, taus, pmin) - outer(taus, taus)
        J <- objects[[1]]$J
        p <- dim(J)[1]
        H <- array(unlist(lapply(objects, function(x) x$Hinv)), c(p, p, m))
        H <- matrix(aperm(H, c(1, 3, 2)), p * m, p) %*% t(chol(J))
        W <- (H %*% t(H)) * (kronecker(Omega, outer(rep(1, p), 
            rep(1, p))))
        coef <- unlist(lapply(objects, function(x) x$coefficients[,1]))
        D <- kronecker(diff(diag(m)), cbind(0,diag(p-1)))
        ndf <- p * (m - 1)
        Tn <- t(D %*% coef) %*% solve(D %*% W %*% t(D), D %*% coef)
        ddf <- n * m - p * (m - 1)
        pvalue <- 1 - pf(Tn, ndf, ddf)
        nobjects <- 1
        tnote1 <- paste("Model: ", models[[1]], "\n", sep = "")
        tnote2 <- paste("Test of Equality of Slopes: tau in { ", 
            paste(taus, collapse = " "), " }\n")
        topnote <- paste(tnote1, tnote2, sep = "")
    }
    table <- data.frame(ndf, ddf, Tn, pvalue)
    dimnames(table)[[2]] <- c("Df", "Resid Df", "F value", "Pr(>F)")
    title <- "Quantile Regression Analysis of Variance Table\n"
    a <- structure(table, heading = c(title, topnote), class = c("anova", 
        "data.frame"))
    print(a)
}
"rq.test.Hajek" <-
function (x0, x1, y, score = "wilcoxon") 
{
    v <- rq(y ~ x0-1, tau = -1)
    r <- ranks(v, score)
    x1hat <- as.matrix(qr.resid(qr(cbind(1, x0)), x1))
    Tn <- as.matrix(t(x1hat) %*% r$ranks)
    Tn <- t(Tn) %*% solve(crossprod(x1hat)) %*% Tn/r$A2
    ndf <- ncol(x1)
    Tn <- Tn/ndf
    ddf <- length(y) - ncol(x0)
    pvalue <- 1 - pf(Tn,ndf,ddf)
    return(Tn,ndf,ddf,pvalue)
}
