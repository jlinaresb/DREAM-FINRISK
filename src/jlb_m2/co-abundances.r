# PIPELINE
# https://biovcnet.github.io/_pages/NetworkScience_igraphcluster.html

# PAPER
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226

# PACKAGE
# https://github.com/zdk123/SpiecEasi

setwd(here::here())
source("requirements.r")
source("src/utils/importPseq.r")
source("src/utils/prepro_functions.r")

# Load data
train <- pseq(subset = "train")
test <- pseq(subset = "test")

# Remove samples
train <- remove_samples(train,
                        remove_nas = TRUE,
                        remove_neg = FALSE)
test <- remove_samples(test,
                       remove_nas = TRUE,
                       remove_neg = FALSE)

# Filter taxa by counts
train <- filter_taxa(train,
                     function(x) sum(x > 2) > (0.8 * length(x)), TRUE)
test <- filter_taxa(test,
                     function(x) sum(x > 2) > (0.8 * length(x)), TRUE)

# Agglomerate by Species
train <- remove_taxa(train, glomby = "Species")
train <- subset_taxa(train, Species != "s__")
test <- remove_taxa(test, glomby = "Species")
train <- subset_taxa(test, Species != "s__")

# Calculate richness
richness_train <- estimate_richness(
                        train,
                        split = TRUE,
                        measures = c("Shannon", "Observed",
                                     "Chao1", "ACE", "Simpson",
                                     "InvSimpson", "Fisher"))
#sample_data(train)$shannon_index <- richness_train$Shannon

richness_test <- estimate_richness(
                        test,
                        split = TRUE,
                        measures = c("Shannon", "Observed",
                                     "Chao1", "ACE", "Simpson",
                                     "InvSimpson", "Fisher"))
#sample_data(test)$shannon_index <- richness_test$Shannon

require(igraph)
require(SpiecEasi)
require(viridis)
require(ggplot2)
require(dplyr)
pargs <- list(seed = 666, ncores = 8, rep.num = 50)
spiec <- spiec.easi(train,
                    method = "glasso",
                    sel.criterion = "stars",
                    lambda.min.ratio = 1e-2,
                    nlambda = 20,
                    pulsar.params = pargs)
cor <- cov2cor(apply(getOptCov(spiec), 2, as.numeric))
weighted_adj_mat <- cor * getRefit(spiec)
grph <- adj2igraph(weighted_adj_mat)

plot(grph, vertex.size = 1,
     vertex.label = NA,
     edge.width = 1,
     layout = layout.circle(grph))

#Remove edges with very low weight 
weight_threshold <- 0.01
grph <- delete.edges(grph, which(abs(E(grph)$weight) < weight_threshold))

#Remove negative edges
grph_pos <- delete.edges(grph, which(E(grph)$weight < 0))
plot(grph_pos,
     vertex.label = NA,
     edge.color = "black",
     layout = layout_with_fr(grph_pos))

#Remove unconnected vertices
grph_pos <- delete.vertices(grph_pos, which(degree(grph_pos) < 1))
plot(grph_pos,
     vertex.label = NA,
     edge.color = "black",
     layout = layout_with_fr(grph_pos))

dd_grph_pos <- degree.distribution(grph_pos)
plot(0:(length(dd_grph_pos) - 1), dd_grph_pos, type = "b",
      ylab = "Frequency", xlab = "Degree", main = "Degree Distributions")


grph_pos_deg <- degree(grph_pos, v = V(grph_pos), mode = "all")
fine <- 500 # this will adjust the resolving power.

#this gives you the colors you want for every point
graphCol <- viridis(fine)[as.numeric(cut(grph_pos_deg, breaks = fine))]

# now plot
plot(grph_pos,
     vertex.color = graphCol,
     edge.color = "black",
     vertex.label = NA,
     layout = layout_with_fr(grph_pos))


grph_pos_bw <- betweenness(grph_pos, directed = FALSE)

#this gives you the colors you want for every point
graphCol <- viridis(fine)[as.numeric(cut(grph_pos_bw, breaks = fine))]

# now plot
plot(grph_pos, vertex.color = graphCol,
     vertex.label = NA,
     edge.color = "black",
     layout = layout_with_fr(grph_pos))


grph_pos_tran <- transitivity(grph_pos, type="local")
grph_pos_tran

#this gives you the colors you want for every point
graphCol <- viridis(fine)[as.numeric(cut(grph_pos_tran, breaks = fine))]

# now plot
plot(grph_pos, vertex.color = graphCol,
     vertex.label = NA,
     edge.color = "black",
     layout = layout_with_fr(grph_pos))


grph_pos_tran_gl <- transitivity(grph_pos, type = "global")

grph_pos_greedy <- cluster_fast_greedy(grph_pos, weights = E(grph_pos)$weight)
modularity(grph_pos_greedy)
sizes(grph_pos_greedy)


colourCount <- length(unique(grph_pos_greedy$membership))
cluster_col <- rainbow(colourCount)[as.numeric(cut(grph_pos_greedy$membership, breaks = colourCount))]

# now plot
plot(grph_pos, vertex.color = cluster_col,
     vertex.label = NA,
     edge.color = "black",
     layout = layout_with_fr(grph_pos))



grph_pos_louvain <- cluster_louvain(grph_pos, weights = E(grph_pos)$weight)
modularity(grph_pos_louvain)
sizes(grph_pos_louvain)

colourCount <- length(unique(grph_pos_louvain$membership))
cluster_col <- rainbow(colourCount)[as.numeric(cut(grph_pos_louvain$membership, breaks = colourCount))]

# now plot
plot(grph_pos, vertex.color = cluster_col,
     vertex.label = NA,
     edge.color = "black",
     layout = layout_with_fr(grph_pos))


V(grph_pos)$cluster <- grph_pos_louvain$membership
vertex_attr(grph_pos, index = V(grph_pos))

ids <- which(sizes(grph_pos_louvain) <= 2)
grph_pos_main_communities <- delete_vertices(grph_pos, which(V(grph_pos)$cluster %in% ids))


nodes <- V(grph_pos_main_communities)$name


taxa <- tax_table(train)
taxa <- as.data.frame(taxa[nodes, ])
taxa$clusters <- V(grph_pos_main_communities)$cluster

View(taxa)
