# PIPELINE
# https://biovcnet.github.io/_pages/NetworkScience_igraphcluster.html

# PAPER
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226

# PACKAGE
# https://github.com/zdk123/SpiecEasi

co_abundances <- function(pseq) {

     pargs <- list(seed = 666, ncores = 8, rep.num = 50)
     spiec <- spiec.easi(pseq,
                         method = "glasso",
                         sel.criterion = "stars",
                         lambda.min.ratio = 1e-2,
                         nlambda = 20,
                         pulsar.params = pargs)

     cor <- cov2cor(apply(getOptCov(spiec), 2, as.numeric))
     weighted_adj_mat <- abs(cor) * getRefit(spiec)
     grph <- adj2igraph(weighted_adj_mat)

     #Remove edges with very low weight 
     weight_threshold <- 0.01
     grph <- delete.edges(grph, which(abs(E(grph)$weight) < weight_threshold))

     #Remove unconnected vertices
     grph <- delete.vertices(grph, which(degree(grph) < 1))

     grph_deg <- degree(grph, v = V(grph), mode = "all")
     fine <- 500 # this will adjust the resolving power.

     #Clustering
     grph_louvain <- cluster_louvain(grph,
                                     weights = E(grph)$weight,
                                     resolution = 4)

     V(grph)$cluster <- grph_louvain$membership
     vertex_attr(grph, index = V(grph))

     nodes <- V(grph)$name

     taxa <- tax_table(pseq)
     taxa <- as.data.frame(taxa[nodes, ])
     taxa$clusters <- V(grph)$cluster

     return(taxa)
}
