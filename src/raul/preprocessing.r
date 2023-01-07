# Load data
print("Importing data ...")
inputdir="./"
train <- pseq(inputdir = inputdir, subset = "train")
test <- pseq(inputdir = inputdir, subset = "test")

####################################################################

# Agglomerate by Species
print("Agglomerating by species ...")
train <- tax_glom(train, taxrank = "Species")
train <- subset_taxa(train, Species != "s__")

test <- tax_glom(test, taxrank = "Species")
test <- subset_taxa(test, Species != "s__")

train_filt <- subset_taxa(train, Phylum %in% c("p__Firmicutes","p__Bacteroidetes"))
test_filt <- subset_taxa(test, Phylum %in% c("p__Firmicutes","p__Bacteroidetes"))
calculate_f_b_ratio <- function(data){
  f_tax <- subset_taxa(data, Phylum == "p__Firmicutes")
  b_tax <- subset_taxa(data, Phylum == "p__Bacteroidetes")
  
  f_abundance <- sample_sums(f_tax)
  b_abundance <- sample_sums(b_tax)
  return(f_abundance / b_abundance)
}

disbiosistrain <- calculate_f_b_ratio(train_filt)
disbiosistest <- calculate_f_b_ratio(test_filt)
sample_data(train)$disbiosis <- disbiosistrain
sample_data(test)$disbiosis <- disbiosistest
########################################################

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

# Caculate GSVA
print("Extract GSVA scores by each cluster ...")
scores <- get_scores(coab_taxa, train, test, method = "gsva")

# Create train and test data
pheno_train <- sample_data(train)
x_train <- data.frame(pheno_train[, -c(8, 9)],
                      richness_train,
                      scores$train)
y_train <- pheno_train[, c("Event_time", "Event")]

pheno_test <- sample_data(test)
x_test <- data.frame(pheno_test[, -c(8, 9)],
                     richness_test,
                     scores$test)
y_test <- pheno_test[, c("Event_time", "Event")]
