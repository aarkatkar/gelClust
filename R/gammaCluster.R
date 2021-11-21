gammaCluster <- function(gm, m, maxiter = 200) {
    A <- -log(gm)
    diag(A) <- 0
    opGamma <- .optimizeGamma(m, A, maxiter)

    return(list(
        labels = opGamma$clusters,
        membership = opGamma$membership
    ))
}


.optimizeGamma <- function(m, A, maxiter) {
    clusters <- Rfast::colMaxs(m, value = FALSE)
    converged <- 0
    for (i in 1:maxiter) {
        d <- eigenMapMatMult(m, A)
        total <- Rfast::rowsums(m) - m
        p <- -pgamma(d, total, lower.tail = FALSE, log.p = TRUE)
        maxCol <- Rfast::colMaxs(p, parallel = TRUE)
        upper <- t(exp(t(p) - maxCol))
        total <- Rfast::colsums(upper, parallel = TRUE)
        m <- t(t(upper) / total)
        clusters0 <- Rfast::colMaxs(m, value = FALSE)
        if (Rfast::all_equals(clusters0, clusters, fast_result = TRUE) &
            i > 10) {
            converged <- 1
            break
        }
        clusters <- clusters0
    }
    return(list(
        clusters = clusters,
        membership = m,
        converged = converged
    ))
}

gammaPlot <- function(gm, clusters) {
    colorList <- scales::hue_pal()(max(clusters))
    u <- as.data.frame(umap::umap(-log1p(-gm), input = "dist")$layout)
    colnames(u) <- c("UMAP1", "UMAP2")
    u$clusters <- as.factor(clusters)
    ggplot2::ggplot() + ggplot2::coord_cartesian() +
        ggplot2::scale_x_continuous() +
        ggplot2::scale_y_continuous() +
        ggplot2::scale_color_hue() +
        ggplot2::layer(
            data = u,
            mapping = ggplot2::aes(x = UMAP1,
                                   y = UMAP2,
                                   color = clusters),
            stat = "identity",
            geom = "point",
            position = ggplot2::position_jitter()
        ) +
        ggplot2::theme_classic()
}

gammaTest <- function(x, gm, rbound, maxiter = 200, plot.show=TRUE, excludeFirst=FALSE) {
    A <- -log(gm)
    diag(A) <- 0

    if(excludeFirst){
        n <- c(3:rbound)
    }else{
        n <- c(2:rbound)
    }
    clustScores <- vapply(n,
                          function(n_)
                              .gammaEval(x, A, n_, maxiter, excludeFirst),
                          numeric(2))
    scores <- clustScores[1,]
    convergence <- clustScores[2,]
    pch <- ifelse(convergence == 0, 4,
                  ifelse(scores == max(scores), 8, 16))
    p.color <- ifelse(convergence == 0,
                      "red",
                      ifelse(scores == max(scores), "blue", "black"))
    if(plot.show){
        plot(
            n,
            scores,
            xlab = "number of clusters",
            ylab = "log-likelihood",
            pch = pch,
            col = p.color
        )
    }
    return(scores)
}

.gammaEval <- function(x, A, n, maxiter, excludeFirst) {
    m <- gelSVD(x, n, excludeFirst)

    opGamma <- .optimizeGamma(m, A, maxiter)
    m <- opGamma$membership
    converged <- opGamma$converged
    mtot <- log(colSums(m^2))
    mpe<- sum(mtot)
    return(c(mpe, converged))
}

subClusters <- function(x, clust, rbound, squared=FALSE, z=NULL, maxiter=200){
    n <- seq_len(dim(clust$membership)[1])
    eg <- eigenGenes(x, clust$membership)
    clust.labels <- vapply(n, function(n_)
        .subCluster(n_,x, clust, eg, rbound, squared, z, maxiter),
        numeric(dim(x)[1]))
    combined <- rowSums(clust.labels)
    mapping <- floor(sort(unique(combined))/rbound)+1
    return(list(
        labels=data.table::frank(combined, ties.method="dense"),
        roots=mapping
    ))
}

.subCluster <- function(n, x, clust, eg, rbound, squared, z, maxiter){
    x.blank <- rep(0, dim(x)[1])
    if(sum(clust$labels==n)==1){
        x.blank[clust$labels==n] <- (n-1)*rbound+1
    }else if(sum(clust$labels==n)==2){
        x.blank[clust$labels==n] <- (n-1)*rbound + c(1,2)
    }else if(sum(clust$labels==n)!=0){
        x.sub <- x[clust$labels==n,]
        gm <- gelMatrix(x.sub, z=rbind(z, eg[n,]), squared=squared)
        scores <- gammaTest(x.sub, gm, rbound, maxiter, plot.show=FALSE, excludeFirst=TRUE)
        n.optimal <- which.max(scores)+2
        me <- gelSVD(x.sub, n.optimal, excludeFirst=TRUE)
        subclusts <- gammaCluster(gm, me, maxiter)
        x.blank[clust$labels==n] <-(n-1)*rbound + subclusts$labels
    }
    return(x.blank)
}

