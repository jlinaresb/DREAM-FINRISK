# Load data
print("Importing data ...")
train <- pseq(inputdir = inputdir, subset = "train")
test <- pseq(inputdir = inputdir, subset = "test")

print(train)
print(test)

# Agglomerate by Species
print("Agglomerating by species ...")
train <- tax_glom(train, taxrank =  "Species")
train <- subset_taxa(train, Species != "s__")
test  <- tax_glom(test, taxrank = "Species")
test  <- subset_taxa(test, Species != "s__")

# Relative abundance
print("Transforming into relative abundance ...")
train_relab <- transform_sample_counts(train, function(x) x / sum(x) )
test_relab  <- transform_sample_counts(test, function(x) x / sum(x) )

# Calculate disbiosis scores
print("Calculating disbiosis scores ...")
train_filt <- subset_taxa(train_relab, Phylum %in% c("p__Firmicutes","p__Bacteroidetes"))
test_filt  <- subset_taxa(test_relab, Phylum %in% c("p__Firmicutes","p__Bacteroidetes"))

calculate_f_b_ratio <- function(data){
  f_tax <- subset_taxa(data, Phylum == "p__Firmicutes")
  b_tax <- subset_taxa(data, Phylum == "p__Bacteroidetes")
  
  f_abundance <- sample_sums(f_tax)
  b_abundance <- sample_sums(b_tax)
  return(f_abundance / b_abundance)
}

disbiosis_train <- calculate_f_b_ratio(train_filt)
disbiosis_test <- calculate_f_b_ratio(test_filt)
sample_data(train)$disbiosis <- disbiosis_train
sample_data(test)$disbiosis <- disbiosis_test

# Relative Abundance of phylos 
print("Calculating relative abundance of each phylo ...")
train_ph <- tax_glom(train_relab, taxrank =  "Phylum")
test_ph <- tax_glom(test_relab, taxrank = "Phylum")

train_ph <- subset_taxa(train_ph, Domain %in% c("k__Archaea", "k__Bacteria"))
test_ph  <- subset_taxa(test_ph, Domain %in% c("k__Archaea", "k__Bacteria"))

train_phylos <- t(train_ph@otu_table@.Data)
colnames(train_phylos) <- tax_table(train_ph)[,2]

test_phylos <- t(test_ph@otu_table@.Data)
colnames(test_phylos) <- tax_table(test_ph)[,2]

# Calculate richness
print("Calculating richness indexes ...")
richness_train <- estimate_richness(
                        train,
                        split = TRUE,
                        measures = c("Shannon", "Observed",
                                     "Chao1", "Simpson",
                                     "InvSimpson", "Fisher"))

richness_test <- estimate_richness(
                        test,
                        split = TRUE,
                        measures = c("Shannon", "Observed",
                                     "Chao1", "Simpson",
                                     "InvSimpson", "Fisher"))

# Filter taxa by counts
print("Filter NZV taxa ...")
train <- filter_taxa(train,
                     function(x) sum(x > 2) > (0.8 * length(x)), TRUE)

# Calculate co-abundances
print("Calculating co-abundances clusters ...")
coab_taxa <- co_abundances(train)
saveRDS(
  coab_taxa,
  file = file.path(
    outputdir, "coabundances.rds"
  ))

# Caculate GSVA
print("Extract GSVA scores by each cluster ...")
scores <- get_scores(coab_taxa, train, test, method = "gsva")

# Create train and test data
pheno_train <- sample_data(train)
x_train <- data.frame(pheno_train[, -c(8, 9)],
                      richness_train,
                      train_phylos,
                      scores$train)
y_train <- pheno_train[, c("Event_time", "Event")]

pheno_test <- sample_data(test)
x_test <- data.frame(pheno_test[, -c(8, 9)],
                     richness_test,
                     test_phylos,
                     scores$test)
y_test <- pheno_test[, c("Event_time", "Event")]


# Check data
# ======
stopifnot(ncol(x_train) == ncol(x_test))
stopifnot(colnames(x_train) == colnames(x_test))
stopifnot(nrow(x_train) == nrow(y_train))
stopifnot(nrow(x_test) == nrow(y_test))


# Creating train and test data
print("Creating train and test data ...")
train <- cbind.data.frame(x_train, y_train)
#test  <- cbind.data.frame(x_test, y_test)
test <- as.data.frame(x_test)


# What do we do with ...
# ======
# NA's values in train and test?
train <- train[complete.cases(train), ]
test <- missRanger(test, pmm.k = 10, seed = 153)
# Negatives survival values in train and test?
train <- train[-which(train$Event_time < 0), ]
train$Event_time <- train$Event_time + 0.1
#test$Event_time[which(test$Event_time < 0)] <- 15

# Removing PrevalentHFAIL (only 0)
train <- train[, -grep("PrevalentHFAIL", colnames(train))]
test <- test[, -grep("PrevalentHFAIL", colnames(test))]


#save(train, test, file = "~/git/DREAM-FINRISK/tmp/data_new.RData")
