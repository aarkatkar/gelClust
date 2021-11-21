eigenGenes <- function(x, m) {
    xm <- apply(x, 1, mean)
    xsd <- apply(x, 1, sd)
    xt <- (x - xm) / xsd
    egenes <- t(apply(m, 1, function(mi) {
          .getEigen(xt, mi)
      }))
}

.getEigen <- function(xt, m) {
    if (sum(m) == 0) {
        egene <- rep(0, dim(xt)[2])
    } else {
        egene <- prcomp(xt * m)$rotation[, 1]
        weightMean <- apply(xt, 2, function(b) {
            weighted.mean(b, m)
        })
        if (cor(egene, weightMean) < 0) {
            egene <- -1 * egene
        }
    }
    return(egene)
}

gelDiff <- function(gm1,
                    gm2,
                    twoTailed = FALSE,
                    unifOut = FALSE) {
    gmDiff <- log(gm2) - log(gm1)
    diag(gmDiff) <- 0
    p <- 0.5 * exp(-abs(gmDiff))
    if (twoTailed) {
        pval <- 2 * p
    } else {
        pval <- ifelse(gmDiff > 0, p, 1 - p)
    }
    N <- dim(gm2)[1]
    if (unifOut) {
        pval <- .forceUniform(pval, lower.tail = TRUE)
    }
    return(pval)
}

.forceUniform <- function(x, lower.tail = FALSE) {
    vals <- x[upper.tri(x)]
    N_ <- length(vals)
    if (lower.tail) {
        rankMat <- data.table::frank(vals) / (N_ + 1)
    } else {
        rankMat <- (N_ + 1 - data.table::frank(vals)) / (N_ + 1)
    }
    p <- diag(0, nrow = dim(x)[1])
    p[upper.tri(p)] <- rankMat
    return(p + t(p))
}
