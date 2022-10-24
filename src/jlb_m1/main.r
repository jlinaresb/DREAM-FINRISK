# TODO:
# - Add guanrank variable to stratify patients
# - Feature selection for each taxon level (PCA, FCBF, filter, wrapper)
# - Trainning models:
#     Option 1): Survival + treatment (like anthra project)
#     Option 2): Binary classification after guanrank stratification
#     Option 3): Survival model


setwd(here::here())
source("requirements.r")
source("src/utils/importPseq.r")
source("src/utils/guanrank.r")

# Import data
train <- pseq(subset = "train")
tt <- train

# Removing artifacts samples with NA survival values
# see https://www.synapse.org/#!Synapse:syn27130803/discussion/threadId=9722
artifacts <- subset_samples(train, PrevalentHFAIL == 1 &
                                   Event == 1 &
                                   Event_time < 0)
nas <- subset_samples(train, is.na(Event_time) | is.na(Event))

to_remove <- c(samples_names(artifacts),
               samples_names(nas))

train <- prune_samples(
            !(sample_names(train) %in% to_remove),
            train)

# Add guanrank score
surv <- data.frame(
    time = get_variable(train, "Event_time"),
    status = get_variable(train, "Event"),
    row.names = rownames(sample_data(train))
)
surv$time <- surv$time + abs(min(surv$time))

grank <- as.data.frame(guanrank(surv))
sample_data(train)$guanrank <- grank$rank


# Select only Bacteria and Archaea kingdom
train <- subset_taxa(train, Domain == "k__Bacteria" |
                            Domain == "k__Archaea")

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
# Log2 transformation
train <- transform_sample_counts(train, function(x) log2(x + 1))

