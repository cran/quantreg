"latex" <- function(x, ...) UseMethod("latex")

"latex.summary.rqs" <-
function (x, transpose = FALSE, caption = "caption goes here.", 
	digits = 3, file = as.character(substitute(x)), ...)
{   
        taus <- function(x) x$tau
        taus <- unlist(lapply(x,taus))
        taus <- format(round(taus, digits))
        coef <- lapply(x,coefficients)
        p <- nrow(coef[[1]])
        k <- ncol(coef[[1]])
        m <- length(taus)
        rlab <- dimnames(coef[[1]])[[1]]
        clab <- taus
        a <- format(round(array(unlist(coef),c(p,k,m)),digits = 3))
        table <- matrix("", p, m)
        for (i in 1:m) {
           for (j in 1:p) {
              if (k == 3) {
                table[j, i] <- paste("$\\underset{(", a[j, 2,
                  i], ",", a[j, 3, i], ")}{", a[j, 1, i], "}$", sep="")
            }
              else if (k == 4) {
                table[j, i] <- paste("$\\underset{(", a[j,2,i] , ")}{",
                        a[j,1, i], "}$",sep="")
               }
            }
        }   
        rowlabel <- "Covariates"
        dimnames(table) <- list(rlab, clab)
        if(transpose) { 
              table <- t(table)
              rowlabel <- "Quantiles"
             }
        latex.table(table, caption = caption, rowlabel = rowlabel, file = file)
        invisible()
}       

"latex.table.rq" <-  
function(x, ...) {
cat("table.rq and related methods are defunct -- see rq, which now accepts a vector of taus.")
}

"plot.table.rq" <- 
function(x, ...){
cat("table.rq() and related methods are defunct --  see ?rq, which now accepts a vector of taus.")
}

"table.rq" <- 
function(x, ...){
cat("table.rq() and related methods are defunct  -- see ?rq, which now accepts a vector of taus.")
}

"plot.summary.rqs" <-
function (x, nrow = 3, ncol = 2, alpha= .1, ols = TRUE, ...) {
# x is an object of class summary.rqs created presumably by  summary.rqs()
# as a list of summary.rq objects the first  task is to unlist stuff:
        taus <- function(x) x$tau
        xx <- unlist(lapply(x,taus))
        coef <- lapply(x,coefficients)
        p <- nrow(coef[[1]])
        k <- ncol(coef[[1]])
        m <- length(xx)
        blab <- dimnames(coef[[1]])[[1]]
        a <- array(unlist(coef),c(p,k,m))
        zalpha <- qnorm(1 - alpha/2)
        par(mfrow = c(nrow,ncol))
        for(i in 1:p){
                if(k == 4){ # standard error methods 
                        b  <- a[i,1,] 
                        bl <- a[i,1,] - a[i,2,]*zalpha
                        bu <- a[i,1,] + a[i,2,]*zalpha
                        }
                else if(k==3){ # rank inversion confidence intervals
                        b  <- a[i,1,] 
                        bl <- a[i,2,] 
                        bu <- a[i,3,] 
                        }
        else stop("summary.rqs components have wrong dimension")
        plot(rep(xx, 2), c(bl,bu), xlab = "", ylab = "", type = "n")
        title(paste(blab[i]),cex = .75)
        polygon(c(xx,rev(xx)),c(bl,rev(bu)),col="LightSkyBlue")
        points(xx, b, cex = .5, pch = "o", col = "blue") 
        lines(xx, b, col = "blue") 
        abline(h=0)
	if(ols){
		abline(h=x$olscoef[i,1],col="red")
		abline(h=x$olscoef[i,1] - x$olscoef[i,2]*zalpha, lty = 2, col="pink")
		abline(h=x$olscoef[i,1] + x$olscoef[i,2]*zalpha, lty = 2, col="pink")
		}
        }
}
"latex.table" <-
function (x, file = as.character(substitute(x)), rowlabel = file, 
    rowlabel.just = "l", cgroup, n.cgroup, rgroup, n.rgroup = NULL, 
    digits, dec, rdec, cdec, append = FALSE, dcolumn = FALSE, cdot = FALSE, 
    longtable = FALSE, table.env = TRUE, lines.page = 40, caption, caption.lot, 
    label = file, double.slash = FALSE, ...) 
{
    nc <- ncol(x)
    nr <- nrow(x)
    if (missing(caption) & !missing(caption.lot)) 
        warning("caption.lot is ignored unless caption is specified")
    if (!longtable & !table.env & !missing(caption)) 
        stop("you must have table.env=TRUE if caption is given")
    if (!missing(digits)) 
        .Options$digits <- digits
    sl <- if (double.slash) 
        "\\\\"
    else "\\"
    rlj <- if (rowlabel.just == "l") 
        "l"
    else "c"
    if (!missing(dec)) {
        if (length(dec) == 1) 
            x <- round(x, dec)
        else {
            if (!is.matrix(dec) || nrow(dec) != nrow(x) || ncol(dec) != 
                ncol(x)) 
                stop("dimensions of dec do not match those of x")
            for (i in 1:nr) for (j in 1:nc) x[i, j] <- round(x[i, 
                j], dec[i, j])
        }
        cx <- format(x)
    }
    else if (!missing(rdec)) {
        cx <- NULL
        for (i in 1:nr) {
            x[i, ] <- round(x[i, ], rdec[i])
            cx <- rbind(cx, format(x[i, ]))
        }
    }
    else if (!missing(cdec)) {
        cx <- NULL
        for (j in 1:nc) {
            x[, j] <- round(x[, j], cdec[j])
            cx <- cbind(cx, format(x[, j]))
        }
    }
    else cx <- format(x)
    cx[is.na(x)] <- ""
    if (dcolumn) 
        sep <- "."
    else {
        #cx <- translate(cx, " ", "~")
        cx <- matrix(chartr(" ", "~", cx), nrow=nr)
        if (cdot) {
            #cx <- translate(cx, "[.]", "\\\\cdot", multichar = TRUE)
            cx <- gsub("[.]", "\\\\cdot", cx)
            cx <- matrix(paste("$", cx, "$", sep = ""), nrow = nr)
            cx[is.na(x)] <- ""
        }
        sep <- "c"
    }
    if (is.null(n.rgroup) && !missing(rgroup)) 
        n.rgroup <- rep(nr/length(rgroup), length(rgroup))
    if (!is.null(n.rgroup) && sum(n.rgroup) != nr) 
        stop("sum of n.rgroup must equal number of rows in x")
    if (!missing(rgroup) && !is.null(n.rgroup) && (length(rgroup) != 
        length(n.rgroup))) 
        stop("lengths of rgroup and n.rgroup must match")
    fi <- paste(file, ".tex", sep = "")
    rowname <- dimnames(x)[[1]]
    if (length(rowname) == 0) {
        rowname <- NULL
        rowlabel <- NULL
        if (!missing(rgroup)) 
            stop("you must have row dimnames to use rgroup")
    }
    #start new file
    if (!append) 
        cat("", file = fi)
    cat("%", deparse(match.call()), "\n%\n", file = fi, append = TRUE)
    if (dcolumn) 
        cat(sl, "newcolumn{.}{D{.}{", sl, "cdot}{-1}}\n", file = fi, 
            append = TRUE)
    if (!is.null(rowlabel)) 
        form <- paste("|", rowlabel.just, "|", sep = "")
    else form <- ""
    f <- paste("|", sep, sep = "", collapse = "")
    if (missing(cgroup)) 
        ff <- c(rep(f, nc), "|")
    else {
        k <- length(cgroup)
        if (missing(n.cgroup)) 
            n.cgroup <- rep(nc/k, k)
        if (sum(n.cgroup) != nc) 
            stop("sum of n.cgroup must equal number of columns")
        if (length(n.cgroup) != length(cgroup)) 
            stop("cgroup and n.cgroup must have same lengths")
        ff <- NULL
        for (i in 1:k) ff <- c(ff, rep(f, n.cgroup[i]), "|")
    }
    form <- paste(form, paste(ff, collapse = ""), sep = "")
    #if(missing(cgroup)) hline <- "" else hline <- paste(sl,"hline",sep="")
    hline <- ""
    if (!missing(caption)) 
        caption <- paste(sl, "caption", if (missing(caption.lot)) 
            NULL
        else paste("[", caption.lot, "]", sep = ""), "{", caption, 
            if (longtable) 
                NULL
            else paste(sl, "label{", label, "}", sep = ""), "}", 
            sep = "")
    if (!longtable) {
        if (table.env) 
            cat(sl, "begin{table}[hptb]\n", sep = "", file = fi, 
                append = TRUE)
        cat(sl, "begin{center}\n", file = fi, sep = "", append = TRUE)
        cat(sl, "begin{tabular}{", form, "} ", sl, "hline", hline, 
            "\n", sep = "", file = fi, append = TRUE)
    }
    else {
        cat(paste(sl, "setlongtables", sep = ""), paste(sl, "begin{longtable}{", 
            form, "}", sep = ""), sep = "\n", file = fi, append = TRUE)
        if (!missing(caption)) 
            cat(caption, sl, sl, "\n", sep = "", file = fi, append = TRUE)
        cat(sl, "hline", hline, "\n", sep = "", file = fi, append = TRUE)
    }
    if (!missing(cgroup)) {
        cgroup <- paste(sl, "bf ", cgroup, sep = "")
        if (is.null(rowlabel)) {
            labs <- c(paste(sl, "multicolumn{", n.cgroup[1], 
                "}{|c||}{", cgroup[1], "}", sep = "", collapse = ""), 
                if (k > 2) paste(sl, "multicolumn{", n.cgroup[c(-1, 
                  -k)], "}{c||}{", cgroup[c(-1, -k)], "}", sep = "") else NULL, 
                paste(sl, "multicolumn{", n.cgroup[k], "}{c|}{", 
                  cgroup[k], "}", sep = ""))
            g <- paste(sl, "hline", sep = "")
        }
        else {
            rowlabel <- paste(sl, "bf ", rowlabel, sep = "")
            labs <- c(paste(sl, "multicolumn{1}{|", rlj, "||}{", 
                rowlabel, "}", sep = ""), paste(sl, "multicolumn{", 
                n.cgroup[-k], "}{c||}{", cgroup[-k], "}", sep = ""), 
                paste(sl, "multicolumn{", n.cgroup[k], "}{c|}{", 
                  cgroup[k], "}", sep = ""))
            g <- paste(sl, "cline{2-", nc + 1, "}", sep = "")
        }
        cat(labs, file = fi, sep = "&", append = TRUE)
        cat(sl, sl, " ", g, "\n", sep = "", file = fi, append = TRUE)
        if (!is.null(rowlabel)) 
            rowlabel <- ""
    }
    collabel <- dimnames(x)[[2]]
    if (is.null(collabel)) 
        collabel <- as.character(1:nc)
    labs <- c(rowlabel, collabel)
    if (missing(cgroup)) {
        if (is.null(rowlabel)) 
            pre <- c(paste(sl, "multicolumn{1}{|c|}{", sep = ""), 
                rep(paste(sl, "multicolumn{1}{c|}{", sep = ""), 
                  nc - 1))
        else pre <- c(paste(sl, "multicolumn{1}{|", rlj, "||}{", 
            sep = ""), rep(paste(sl, "multicolumn{1}{c|}{", sep = ""), 
            nc))
    }
    else {
        if (is.null(rowlabel)) {
            pre <- NULL
            j <- 0
            for (i in 1:k) {
                if (n.cgroup[i] > 1) {
                  g <- rep(paste(sl, "multicolumn{1}{c|}{", sep = ""), 
                    n.cgroup[i] - 1)
                  if (j == 0) 
                    g[1] <- paste(sl, "multicolumn{1}{|c|}{", 
                      sep = "")
                  pre <- c(pre, g)
                }
                j <- j + n.cgroup[i]
                if (j == 1) 
                  g <- paste(sl, "multicolumn{1}{|c||}{", sep = "")
                else if (j < nc) 
                  g <- paste(sl, "multicolumn{1}{c||}{", sep = "")
                else g <- paste(sl, "multicolumn{1}{c|}{", sep = "")
                pre <- c(pre, g)
            }
        }
        else {
            pre <- paste(sl, "multicolumn{1}{|", rlj, "||}{", 
                sep = "")
            j <- 0
            for (i in 1:k) {
                pre <- c(pre, rep(paste(sl, "multicolumn{1}{c|}{", 
                  sep = ""), n.cgroup[i] - 1))
                j <- j + n.cgroup[i]
                if (j < nc) 
                  g <- paste(sl, "multicolumn{1}{c||}{", sep = "")
                else g <- paste(sl, "multicolumn{1}{c|}{", sep = "")
                pre <- c(pre, g)
            }
        }
    }
    labs <- paste(pre, labs, "}", sep = "")
    cat(labs, file = fi, sep = "&", append = TRUE)
    cat(sl, sl, " ", sl, "hline", hline, "\n", sep = "", file = fi, 
        append = TRUE)
    if (longtable) {
        if (missing(caption)) 
            cat(sl, "endhead\n", sl, "hline", sl, "endfoot\n", 
                sep = "", file = fi, append = TRUE)
        else {
            cat(sl, "endfirsthead\n", sep = "", file = fi, append = TRUE)
            if (!missing(caption)) 
                cat(sl, "caption[]{\\em (continued)} ", sl, sl, 
                  "\n", sep = "", file = fi, append = TRUE)
            cat(sl, "hline", hline, "\n", sep = "", file = fi, 
                append = TRUE)
            cat(labs, file = fi, sep = "&", append = TRUE)
            cat(sl, sl, " ", sl, "hline", hline, "\n", sl, "endhead", 
                sl, "hline", sl, "endfoot\n", sep = "", file = fi, 
                append = TRUE)
            cat(sl, "label{", label, "}\n", sep = "", file = fi, 
                append = TRUE)
        }
    }
    if (is.null(n.rgroup)) 
        rg.end <- 0
    else {
        rg.end <- cumsum(n.rgroup)
        rg.start <- rg.end - n.rgroup + 1
        if (missing(rgroup)) 
            rgroup <- rep("", length(n.rgroup))
        else rgroup <- paste("{", sl, "bf ", rgroup, "}", sep = "")
    }
    linecnt <- 0
    for (i in 1:nr) {
        if (!missing(rgroup)) {
            k <- rg.start == i
            if (any(k)) {
                j <- (1:length(n.rgroup))[k]
                if (longtable && linecnt > 0 && (linecnt + n.rgroup[j] + 
                  (n.rgroup[j] > 1)) > lines.page) {
                  cat(sl, "newpage\n", sep = "", file = fi, append = TRUE)
                  linecnt <- 0
                }
                if (n.rgroup[j] > 1) {
                  cat(rgroup[j], rep("", nc), file = fi, sep = "&", 
                    append = TRUE)
                  linecnt <- linecnt + 1
                  cat(sl, sl, "\n", sep = "", file = fi, append = TRUE)
                }
                l <- rg.start[j]:rg.end[j]
                if (length(l) > 1) 
                  rowname[l] <- paste("~~", rowname[l], sep = "")
                else rowname[l] <- paste("{", sl, "bf ", rowname[l], 
                  "}", sep = "")
            }
        }
        else if (longtable && linecnt > 0 && (linecnt + 1 > lines.page)) {
            cat(sl, "newpage\n", sep = "", file = fi, append = TRUE)
            linecnt <- 0
        }
        cat(c(rowname[i], cx[i, ]), file = fi, sep = "&", append = TRUE)
        linecnt <- linecnt + 1
        if (i < nr && any(rg.end == i)) 
            g <- paste(sl, "hline", sep = "")
        else g <- ""
        cat(sl, sl, " ", g, "\n", sep = "", file = fi, append = TRUE)
    }
    cat(sl, "hline", hline, "\n", sep = "", file = fi, append = TRUE)
    if (longtable) 
        cat(sl, "end{longtable}\n", sep = "", file = fi, append = TRUE)
    else {
        cat(sl, "end{tabular}\n", sep = "", file = fi, append = TRUE)
        if (!missing(caption)) 
            cat(sl, "vspace{3mm}\n", sep = "", file = fi, append = TRUE)
        cat(caption, "\n", file = fi, append = TRUE)
        cat(sl, "end{center}\n", sep = "", file = fi, append = TRUE)
        if (table.env) 
            cat(sl, "end{table}\n", sep = "", file = fi, append = TRUE)
    }
    invisible()
}
"table.rq" <-
function (x, ...) 
stop("table.rq now defunct, rq() now accepts vector tau argument.  See ?rq.")
