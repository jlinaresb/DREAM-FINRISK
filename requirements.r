if (!("phyloseq" %in% installed.packages())){BiocManager::install("phyloseq",ask = F,update = F)}
if (!("mixOmics" %in% installed.packages())){BiocManager::install("mixOmics",ask = F,update = F)}
if (!("survcomp" %in% installed.packages())){BiocManager::install("survcomp",ask = F,update = F)}
if (!("biospear" %in% installed.packages())){install.packages("biospear")}

require(phyloseq)
require(biospear)
