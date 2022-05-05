gelMatrix <- function(x,
                      phen = NULL,
                      z = NULL,
                      squared = FALSE,
                      partial = FALSE,
                      output = "unif") {
    if (output != "unif" & output != "orig" & output != "cor") {
        warning('output must be "unif", "orig", or "cor"')
        return(NA)
    }

    # Get correlations
    correlations <- .getCorrel(x, phen, z, partial)
    corMat <- correlations$corMat
    if (output == "cor") {
        return(corMat)
    }

    n_phen <- correlations$n_phen
    n_cov <- correlations$n_cov
    # Get p-values
    if (partial | output == "unif") {
        p <- .unifPvalue(corMat, n_phen, dim(x)[2], squared)
    } else {
        p <- .tdistPvalue(corMat, n_phen, n_cov, dim(x)[2], squared)
    }
    return(p)
}

.unifPvalue <- function(corMat, n_phen, n_samp, squared) {
    N <- dim(corMat)[1]
    if (n_phen == 0) {
        rxy <- corMat
    } else {
        rxy <- corMat[(n_phen + 1):N, (n_phen + 1):N]
        dm <- corMat[1:n_phen, 1:N]
        suppressWarnings(tapp <- dm * sqrt((n_samp - 2) / (1 - dm^2)))
        if (squared) {
            papp <- 2 * pt(-abs(tapp), n_samp - 2)
        } else {
            papp <- pt(-tapp, n_samp - 2)
        }
    }
    if (squared) {
        rxy <- rxy^2
    }
    p <- .forceUniform(rxy)
    if (n_phen != 0) {
        p <- rbind(papp[1:n_phen, (n_phen + 1):N], p)
        p <- cbind(t(papp), p)
    }
    return(p)
}

.tdistPvalue <- function(corMat, n_phen, n_cov, n_samp, squared) {
    N <- dim(corMat)[1]
    if (n_phen == 0) {
        rxy <- corMat
    } else {
        rxy <- corMat[(n_phen + 1):N, (n_phen + 1):N]
        dm <- corMat[1:n_phen, 1:N]
        suppressWarnings(tapp <- dm * sqrt((n_samp - 2) / (1 - dm^2)))
        if (squared) {
            papp <- 2 * pt(-abs(tapp), n_samp - 2)
        } else {
            papp <- pt(-tapp, n_samp - 2)
        }
    }
    rxy <- ifelse(rxy >= 1, 1, rxy)
    suppressWarnings(t <-
        rxy * sqrt((n_samp - 2 - n_cov) / (1 - rxy^2)))
    if (squared) {
        p <- 2 * pt(-abs(t), n_samp - 2 - n_cov)
    } else {
        p <- pt(-t, n_samp - 2 - n_cov)
    }
    if (n_phen != 0) {
        p <- rbind(papp[1:n_phen, (n_phen + 1):N], p)
        p <- cbind(t(papp), p)
    }
    p <- matrix(c(p), nrow = N)
    return(p)
}


.getCorrel <- function(x, phen, z, partial) {
    if (partial) {
        corMat <- corpcor::pcor.shrink(t(x), verbose = FALSE)
        n_cov <- 0
    } else if (is.null(z)) {
        corMat <- coop::pcor(t(x), use = "complete.obs")
        n_cov <- 0
    } else {
        if (is.null(dim(z))) {
            z <- t(z)
        }
        corMat <- .covCor(x, z)
        n_cov <- dim(z)[1]
    }
    if (!is.null(phen)) {
        if (is.null(dim(phen))) {
            phen <- t(phen)
        }
        corMat <- .addPhen(corMat, x, phen)
        n_phen <- dim(phen)[1]
    } else {
        n_phen <- 0
    }
    return(list(
        corMat = corMat,
        n_phen = n_phen,
        n_cov = n_cov
    ))
}

.covCor <- function(x, z) {
    N <- dim(x)[1]
    xtemp <- rbind(x, z)
    rxy <- coop::pcor(t(xtemp), use = "complete.obs")
    for (i in 1:dim(z)[1]) {
        rxz <- rxy[N + i, ]
        suppressWarnings(sqrxz <- ifelse(rxz > 1, 0, sqrt(1 - rxz^2)))
        denom <- sqrxz %*% t(sqrxz)
        rxy <- ifelse(denom == 0, 1, (rxy - rxz %*% t(rxz)) / denom)
    }
    rxy <- rxy[1:N, 1:N]
    return(rxy)
}

.addPhen <- function(rxy, x, phen) {
    N <- dim(x)[1]
    nphen <- dim(phen)[1]
    dm <- cor(t(phen), t(rbind(phen, x)), use = "complete.obs")
    rxy <- rbind(dm[1:nphen, (nphen + 1):(N + nphen)], rxy)
    rxy <- cbind(t(dm), rxy)
    return(rxy)
}

gelSVD <- function(x, k, excludeFirst=FALSE, unsigned=FALSE){
    x.standard <- (x - rowMeans(x))/apply(x, 1, sd)
    sv <- svd(x.standard, nu=k, nv=k)
    if (unsigned){
        vexp <- sv$d[1:k]^2
        V <- t(sv$v)
    }else{
        ppos <- colSums(ifelse(sv$u<0,0,sv$u))/colSums(abs(sv$u))
        vexp <- c(ppos*(sv$d[1:k]^2), (1-ppos)*(sv$d[1:k]^2))
        V <- rbind(t(sv$v), -t(sv$v))
    }

    rexp <- rank(vexp, ties.method="first")
    if (excludeFirst & unsigned){
        eg <- V[which.max(vexp),]
        selected <- which(rexp<k)
    }else if (excludeFirst & !unsigned){
        eg <- V[which.max(vexp),]
        selected <- which(rexp>k&rexp<(2*k))
    }
    else if (unsigned){
        eg <- NULL
        selected <- which(rexp>0)
    }else{
        eg <- NULL
        selected <- which(rexp>k)
    }
    V <- V[selected,]
    k <- length(selected)
    p <- -log(gelMatrix(as.matrix(x.standard), phen=V)[1:k,(k+1):(dim(x)[1]+k)])
    maxCol <- Rfast::colMaxs(p, parallel = TRUE)
    upper <- t(exp(t(p) - maxCol))
    total <- Rfast::colsums(upper, parallel = TRUE)
    m <- t(t(upper) / total)
    return(m)
}
