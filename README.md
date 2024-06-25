
<!-- README.md is generated from README.Rmd. Please edit that file -->

This is a repository for the manuscript Tsang et al. (Submitted). Some
datasets are private thus were not uploaded here. Please refer to
supplementary table S3 of the manuscript and contact the data owners for
the data.

Code required to analyze the diversity data has been uploaded, but note
that the results wont be the same with the manuscript as some private
datasets were not released.

Three data files are available: cleaned_alpha.csv, cleaned_beta.csv, and
cleaned_genus.csv

## Code

In the code folder, there are 5 scripts. Please run the scripts in
numerical order

| Script | Description                                                                      |
|--------|----------------------------------------------------------------------------------|
| Obj1a  | Alpha diversity analyses                                                         |
| Obj1b  | Beta diversity analyses                                                          |
| Obj2   | Analyzing the scaling relationships between taxonomic and phylogenetic diversity |
| Obj4   | Genus analyses                                                                   |
| sim    | Simulations about mixed models (Text S2)                                         |

## Alpha diversity analyses

Please use the cleaned_alpha.csv for these analyses.

| Variable               | Description                                                |
|------------------------|------------------------------------------------------------|
| studyID                | Study identity                                             |
| Lat                    | Latitude                                                   |
| Long                   | Longitude                                                  |
| Exclude.AM.Est.Shannon | Estimated Shannon diversity after excluding Apis mellifera |
| qPD.S                  | Estimated phylogenetic diversity when q (Hill number) = 1  |
| Agr1                   | Agricultural site = 1, non-agricultural site = 0           |
| Urb1                   | Urban site = 1, non-urban site = 0                         |
| abundanceMethod        | Sampling methods                                           |
| PC1                    | Climatic PC1 score                                         |
| PC2                    | Climatic PC2 score                                         |
| SC                     | Sampling completeness                                      |

## Beta diversity analyses

Please use the cleaned_beta.csv for these analyses.

| Variable          | Description                                                          |
|-------------------|----------------------------------------------------------------------|
| study             | Study identity                                                       |
| T_beta            | Taxonomic beta diversity                                             |
| P_beta            | Phylogenetic beta diversity                                          |
| metrics           | Btotal (Total beta), Bab (abundance difference), or Bturn (turnover) |
| LC                | Land use (Natural/Agricultural/Urban)                                |
| mean.spatial.dist | Mean spatial distance                                                |
| SC.diff.mean      | Mean sampling completeness differences between assemblages           |
| abundanceMethod   | Sampling methods                                                     |
| mean_PC1          | Mean climatic PC1 score                                              |
| mean_PC2          | Mean climatic PC2 score                                              |
| site              | The number of site in the study                                      |

## genus analyses

Please use the genus.csv for these analyses.

| Variable         | Description                                      |
|------------------|--------------------------------------------------|
| studyID          | Study identity                                   |
| genus            | Genus identity                                   |
| Lat              | Latitude                                         |
| Long             | Longitude                                        |
| abundance        | Abundance value                                  |
| Agriculture      | Agricultural site = 1, non-agricultural site = 0 |
| Urban            | Urban site = 1, non-urban site = 0               |
| abundanceMethod  | Sampling Method                                  |
| SamplingEffort   | Sampling effort at the site                      |
| PC1              | Climatic PC1 score                               |
| PC2              | Climatic PC2 score                               |
| LC.author.coarse | Land use change (single column                   |
