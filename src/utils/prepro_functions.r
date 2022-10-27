remove_samples <- function(pseq,
                           remove_nas = TRUE,
                           remove_neg = TRUE) {
    # Removing artifacts samples with NA survival values and negative times
    # see https://www.synapse.org/#!Synapse:syn27130803/discussion/threadId=9722
    artifacts <- subset_samples(pseq, PrevalentHFAIL == 1 &
                                      Event == 1 &
                                      Event_time < 0)
    nas <- subset_samples(pseq, is.na(Event_time) | is.na(Event))
    negatives <- subset_samples(pseq, Event_time < 0)

    res <- prune_samples(
                !(sample_names(pseq) %in% sample_names(artifacts)),
                pseq)
    # Remove NAs?
    if (remove_nas == TRUE) {
        res <- prune_samples(
            !(sample_names(res) %in% sample_names(nas)),
        res)
        print("Samples with NA's in time or event have been removed!")
    }
    # Remove negatives?
    if (remove_neg == TRUE) {
        res <- prune_samples(
            !(sample_names(res) %in% sample_names(negatives)),
        res)
        print("Samples with times negatives have been removed!")
    }

    return(res)
}

remove_taxa <- function(pseq, glomby = "Species") {
    # Select only Bacteria and Archaea kingdom
    #pseq <- subset_taxa(pseq, Domain == "k__Bacteria" |
    #                         Domain == "k__Archaea")
    # Remove taxa not seen more than 3 times in at least 20% of the samples
    res <- filter_taxa(pseq,
                        function(x) sum(x > 3) > (0.2 * length(x)), TRUE)
    # Agglomerate by species
    res <- tax_glom(res, glomby)
    return(res)
}

norm <- function(pseq, filter = TRUE) {
    # Standardize abundances to the median sequencing depth
    total <- median(sample_sums(pseq))
    standf <- function(x, t = total) round(t * (x / sum(x)))
    res <- transform_sample_counts(pseq,
                                   standf)
    if (filter == TRUE) {
        # Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
        res <- filter_taxa(res,
                           function(x) sd(x) / mean(x) > 3.0, TRUE)
    }
    # Log2 transformation
    res <- transform_sample_counts(res, function(x) log2(x + 1))
    return(res)
}