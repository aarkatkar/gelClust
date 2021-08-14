gammaCluster <- function(gm, n, maxiter=200){
  initSeeds <- .initialSeeds(gm, n)
  m <- initSeeds$membership
  A <- initSeeds$A
  D <- initSeeds$D
  m <- .initialMembership(m, D, n)
  opGamma <- .optimizeGamma(m, A, maxiter)

  return(list(labels=opGamma$clusters, membership=opGamma$membership))
}

.initialSeeds <- function(gm, n){
  A <- -log(gm)
  diag(A) <- 0
  D <- -log1p(-gm)
  diag(D) <- 0
  N <- dim(A)[1]
  bld <- cluster::pam(D, n, diss=TRUE, do.swap=FALSE)
  selected <- bld$medoids

  m <- matrix(rep(0, N*n), nrow=n)
  for (i in seq_len(n)){
    m[i,selected[i]]=1
  }
  return(list(membership=m, A=A, D=D))
}

.optimizeGamma <- function(m, A, maxiter){
  clusters <- Rfast::colMaxs(m, value=FALSE)
  converged <- 0
  for (i in 1:maxiter){
    d <- eigenMapMatMult(m, A)
    total <- Rfast::rowsums(m)-m
    p <- -pgamma(d, total, lower.tail=FALSE, log.p=TRUE)
    maxCol <- Rfast::colMaxs(p, parallel=TRUE)
    upper <- t(exp(t(p)-maxCol))
    total <- Rfast::colsums(upper, parallel=TRUE)
    m <- t(t(upper)/total)
    clusters0 <- Rfast::colMaxs(m, value=FALSE)
    if (Rfast::all_equals(clusters0, clusters, fast_result=TRUE) & i > 10){
      converged <- 1
      break
    }
    clusters <- clusters0
  }
  return(list(clusters=clusters, membership=m, converged=converged))
}

.initialMembership <- function(m, D, n){
  d <- -n*eigenMapMatMult(m, D)
  upper <- exp(d-max(d))
  total <- rowSums(t(upper))
  m <- t(t(upper)/total)
  return(m)
}

gammaPlot <- function(gm, clusters){
  colorList <- scales::hue_pal()(max(clusters))
  u <- as.data.frame(umap::umap(-log1p(-gm), input="dist")$layout)
  colnames(u) <- c("UMAP1", "UMAP2")
  u$clusters <- as.factor(clusters)
  ggplot2::ggplot() + ggplot2::coord_cartesian() +
    ggplot2::scale_x_continuous() +
    ggplot2::scale_y_continuous() +
    ggplot2::scale_color_hue() +
    ggplot2::layer(data=u,
                   mapping=ggplot2::aes(x=UMAP1, y=UMAP2, color=clusters),
                   stat="identity",
                   geom="point",
                   position=ggplot2::position_jitter()
    )
}

gammaTest <- function(start, end, gm, maxiter=200){
  initSeeds <- .initialSeeds(gm, end)
  m <- initSeeds$membership
  A <- initSeeds$A
  D <- initSeeds$D

  n <- c(start:end)
  clustScores <- vapply(n,
                        function(n_) .gammaEval(m, A, D, n_, maxiter),
                        numeric(2))
  scores <- clustScores[1,]
  convergence <- clustScores[2,]
  pch <- ifelse(convergence==0,4,
                ifelse(scores==min(scores), 8, 16))
  p.color <- ifelse(convergence==0,"red",
                    ifelse(scores==min(scores), "blue", "black"))
  plot(n, scores, xlab="n", ylab="Modified Partition Entropy",
       pch=pch, col=p.color)

}

.gammaEval <- function(m, A, D, n, maxiter){
  m <- .initialMembership(m[seq_len(n),], D, n)
  opGamma <- .optimizeGamma(m, A, maxiter)
  m <- opGamma$membership
  converged <- opGamma$converged
  N <- dim(A)[1]
  partEnt <- -1/N*sum(ifelse(m==0, 0, m*log(m)))
  mpe <- N*partEnt/(N-n)
  return(c(mpe, converged))
}

