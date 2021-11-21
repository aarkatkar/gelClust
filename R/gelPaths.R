expMST <- function(gm, cutoff = 0.05, min.size = 3) {
    g <- igraph::graph.adjacency(-log1p(-gm), weighted = TRUE)
    mst <- igraph::mst(g)
    drop <- which(igraph::E(mst)$weight > cutoff)
    filtered <- igraph::delete_edges(mst, igraph::E(mst)[drop])
    clust <- igraph::components(filtered)
    clust.labels <- clust$membership
    clust.size <- clust$csize
    disconnected <- which(clust.size[clust.labels] < min.size)
    clust.labels[disconnected] <- 0
    clusterOld <- c(0, unique(clust.labels[clust.labels != 0]))
    clust.labels <- match(clust.labels, clusterOld) - 1
    return(clust.labels)
}

minPath <- function(s, t, gm, clusters,
                    r = 1,
                    plotpath = TRUE,
                    nodeNames = NULL) {
    if (is.null(nodeNames)) {
        nodeNames <- seq_len(dim(gm)[1])
    }
    s <- which(nodeNames == s)
    t <- which(nodeNames == t)

    n_phen <- dim(gm)[1] - length(clusters)
    odist <- -log1p(-gm)
    tempgm <- odist - log(r)
    diag(tempgm) <- 0
    g <- igraph::graph_from_adjacency_matrix(tempgm, weighted = TRUE)
    sp <- igraph::shortest_paths(g, s, t, output = "both")

    pvals <- igraph::E(g)[sp$epath[[1]]]$weight
    sp_vert <- sp$vpath[[1]]
    nodeNames <- nodeNames[sp_vert]
    nodeClusters <- c(rep(0, n_phen), clusters)[sp_vert] + 1
    gdf <- NULL
    if (plotpath) {
        .plotPath(sp_vert, clusters, nodeNames, nodeClusters)
    }
    return(list(
        pvals = pvals, clusters = nodeClusters - 1,
        names = nodeNames
    ))
}

.plotPath <- function(sp_vert, clusters, nodeNames, nodeClusters) {
    gdf <- data.frame(from = seq_len(length(sp_vert) - 1),
                      to = c(2:length(sp_vert)))
    gdf <- igraph::graph.data.frame(gdf)


    hpal <- c("white", scales::hue_pal()(max(clusters)))
    spal <- c("none", rep("circle", max(clusters)))


    plot(gdf,
        layout = igraph::layout_as_tree, vertex.label = nodeNames,
        vertex.color = hpal[nodeClusters],
        edge.arrow.size = 0.5, vertex.shape = spal[nodeClusters]
    )
}

distBoot <- function(x, partial = FALSE, squared = TRUE, output = "orig") {
    n <- dim(x)[2]
    xShuf <- t(apply(x, 1, sample))
    gmShuf <- gelMatrix(xShuf,
                        squared = squared, partial = partial, output = output)
    D <- -log1p(-gmShuf)
    diag(D) <- 0

    d <- Rfast::floyd(D)
    du <- d[upper.tri(d)]
    return(du)
}

pBoot <- function(d, z) {
    delta <- z - d
    p <- pgamma(delta[delta > 0], 2, lower.tail = FALSE, log.p = TRUE)
    log.p <- sum(p)
    return(1 - exp(log.p))
}
