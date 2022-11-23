if (!("phyloseq" %in% installed.packages())){BiocManager::install("phyloseq",ask = F,update = F)}
if (!("mixOmics" %in% installed.packages())){BiocManager::install("mixOmics",ask = F,update = F)}
if (!("survcomp" %in% installed.packages())){BiocManager::install("survcomp",ask = F,update = F)}
if (!("biospear" %in% installed.packages())){install.packages("biospear")}
if (!("devtools" %in% installed.packages())){install.packages("devtools")}
if (!("igraph" %in% installed.packages())){install.packages("igraph")}
if (!("SpiecEasi" %in% installed.packages())){install_github("zdk123/SpiecEasi")}
if (!("dplyr" %in% installed.packages())){install.packages("dplyr")}
if (!("missRanger" %in% installed.packages())){install.packages("missRanger")}
if (!("GSVA" %in% installed.packages())){BiocManager::install("GSVA",ask = F,update = F)}

require(phyloseq)
require(biospear)
require(igraph)
require(devtools)
require(SpiecEasi)
require(dplyr)
require(GSVA)
require(missRanger)
