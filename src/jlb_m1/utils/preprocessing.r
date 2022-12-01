# Load data
print("Importing data ...")
train <- pseq(inputdir = inputdir, subset = "train")
test <- pseq(inputdir = inputdir, subset = "test")

# Agglomerate by Species
print("Agglomerating by species ...")
train <- tax_glom(train, taxrank =  "Species")
train <- subset_taxa(train, Species != "s__")
test <- tax_glom(test, taxrank = "Species")
test <- subset_taxa(test, Species != "s__")

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
