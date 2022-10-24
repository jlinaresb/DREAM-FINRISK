relative_abundance <- function(x) {
                        return(x / sum(x))
                        }

median_sequencing_depth <- function(x, t) {
                            return(round(t * (x / sum(x))))
                            }