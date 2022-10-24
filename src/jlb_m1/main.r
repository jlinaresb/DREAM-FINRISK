setwd(here::here())
source("requirements.r")
source("src/utils/importPseq.r")
source("src/utils/normalizations.r")

# Import data
train <- pseq(subset = "train")

# Removing artifacts
# see https://www.synapse.org/#!Synapse:syn27130803/discussion/threadId=9722
artifacts <- subset_samples(train, PrevalentHFAIL == 1 &
                                   Event == 1 &
                                   Event_time < 0)
train <- prune_samples(
            !(sample_names(train) %in% sample_names(artifacts)),
            train)

# Remove taxa not seen more than 3 times in at least 20% of the samples
train <- filter_taxa(train,
                     function(x) sum(x > 3) > (0.2 * length(x)), TRUE)
# Standardize abundances to the median sequencing depth
total <- median(sample_sums(train))
standf <- function(x, t = total) round(t * (x / sum(x)))
train <- transform_sample_counts(train,
                                 standf)
# Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
train <- filter_taxa(train,
                     function(x) sd(x) / mean(x) > 3.0, TRUE)
