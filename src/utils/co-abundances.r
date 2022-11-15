# PIPELINE
# https://biovcnet.github.io/_pages/NetworkScience_igraphcluster.html

# PAPER
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226

# PACKAGE
# https://github.com/zdk123/SpiecEasi

co_abundances <- function(train, test, method = "gsva") {

     require(igraph)
     require(SpiecEasi)
     require(viridis)
     require(ggplot2)
     require(dplyr)
     require(GSVA)

     pargs <- list(seed = 666, ncores = 8, rep.num = 50)
     spiec <- spiec.easi(train,
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

     #Degree distributions
     dd_grph <- degree.distribution(grph)
     plot(0:(length(dd_grph) - 1), dd_grph, type = "b",
          ylab = "Frequency", xlab = "Degree", main = "Degree Distributions")


     grph_deg <- degree(grph, v = V(grph), mode = "all")
     fine <- 500 # this will adjust the resolving power.

     #this gives you the colors you want for every point
     graphCol <- viridis(fine)[as.numeric(cut(grph_deg, breaks = fine))]

     # now plot
     plot(grph,
          vertex.color = graphCol,
          edge.color = "black",
          vertex.label = NA,
          layout = layout_with_fr(grph))


     grph_bw <- betweenness(grph, directed = FALSE)

     #Clustering
     grph_louvain <- cluster_louvain(grph, weights = E(grph)$weight, resolution = 4)
     print(paste0("Modularity:", modularity(grph_louvain)))
     sizes(grph_louvain)

     colourCount <- length(unique(grph_louvain$membership))
     cluster_col <- rainbow(colourCount)[as.numeric(cut(grph_louvain$membership, breaks = colourCount))]

     # now plot
     plot(grph, vertex.color = cluster_col,
          vertex.label = NA,
          edge.color = "black",
          layout = layout_with_fr(grph))


     V(grph)$cluster <- grph_louvain$membership
     vertex_attr(grph, index = V(grph))

     nodes <- V(grph)$name

     taxa <- tax_table(train)
     taxa <- as.data.frame(taxa[nodes,])
     taxa$clusters <- V(grph)$cluster

     #TRAIN
     tt <- otu_table(train)@.Data
     tt <- tt[match(rownames(taxa), rownames(tt)), ]
     tt <- t(apply(tt, 1, function(x) log2(x + 1)))

     kk <- unique(taxa$clusters)
     signatures <- list()
     for (i in seq_along(kk)){
          signatures[[i]] <- rownames(taxa[which(taxa$clusters == kk[i]), ])
     }
     names(signatures) <- paste0("cluster_", kk)

     gsva_train <- gsva(tt,
               signatures,
               method = method,
               kcdf = "Gaussian")
     gsva_train <- as.data.frame(t(gsva_train))

     #TEST
     tt <- otu_table(test)@.Data
     tt <- tt[match(rownames(taxa), rownames(tt)), ]
     tt <- t(apply(tt, 1, function(x) log2(x + 1)))

     gsva_test <- gsva(tt,
               signatures,
               method = method,
               kcdf = "Gaussian")
     gsva_test <- as.data.frame(t(gsva_test))


     res <- list(train = gsva_train,
                 test = gsva_test)

     return(res)
}