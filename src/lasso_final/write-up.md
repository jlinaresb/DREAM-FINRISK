# Prediction of Heart Failure risk using biomarker selection in penalized regression models from metagenomic data
Jose Liñares-Blanco <sup>1</sup>, Samuel Pérez-Fernández <sup>1</sup>, Marina Vargas <sup>1</sup>, Ivan Ellson-Lancho <sup>1</sup>, Juan A. Villatoro <sup>1</sup>, Raul López-Domínguez <sup>1</sup>, Jordi Martorell-Marugan <sup>1</sup>, Daniel Toro-Domínguez <sup>1</sup>, Adrián García-Moreno <sup>1</sup> & Pedro Carmona <sup>1</sup>. 

<sup>1</sup>  *GENYO. Centre for Genomics and Oncological Research: Pfizer, University of Granada, Andalusian Regional Goverment, PTS Granada, Avenida de la Ilustración 114, 18016, Granada, Spain*

## Methods

### Data preprocessing
Both train and test data were imported as phyloseq [[1]](https://joey711.github.io/phyloseq/) R object, which includes OTU, taxonomic and metadata table. 

In order to reduce the high dimensionality from metagenomic features, we extracted several characteristics of each patient. First, we agglomerate OTUs to species level and those undefined species were removed (i.e. when species name was *s__*). Second, we calculated several richness scores from agglomerated species, such as Shannon, Obseved, Chao1, Simpson, InvSimpson and Fisher index. In addition, we used the ratio of the relatives abundances between Firmicutes and Bacteroidete phylum, as an indicator of disbiosis. In addition, those Phylos belonging to the kingdoms Archaea and Bacteria were agglomerated and the relative abundance of each was calculated. 


### Co-ocurrence network analysis
To examine co-ocurrence patterns in the metagenomic data, we grouped all taxa at species level by summatory and looked for significant correlations between different species sequence abundances using SparCC method (from the [(SpiecEasi)](https://github.com/zdk123/SpiecEasi)  package in R). Only those genera with significant correlations (both positives and negatives, with p < 0.01) by bootstrapping 50 times were included in the networks. We created a network from all patients in train dataset, where a node represents a specie, and an edge between two species represents a significant correlation in abundances between the pair. To calculate p-values for the network statistics, we performed 500 permutations of the network. After cleaning the network, we aim to finding community structure between species. We ran an unsupervised two-step method for community detection called [Louvain](https://arxiv.org/abs/0803.0476), setting 4 as resolution value. 


### Microbial set variation analysis
Using the [GSVA]() R package, we performed a pathway-centric analysis of metagenomic data. Each cluster obtained from co-ocurrence analyses was reduced to a single vector, used as a predictor in our survival model.


### Biomarker selection
To perform variable selection we used [biospear](https://cran.r-project.org/web/packages/biospear/index.html) R package that implements a large panel of approaches to develop and evaluate a prediction model within a high-dimensional Cox regression setting, and to estimate expected survival at a given time point. Specifically, we selected variables via lasso model. 

All microbiome features were used as biomarkers.

## References
