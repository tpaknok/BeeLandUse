
<!-- README.md is generated from README.Rmd. Please edit that file -->

This is a repository for the manuscript Tsang et al. (Submitted). Some
datasets are private thus were not uploaded here (and therefore the numbers will be a bit different from the main text). Please refer to
supplementary table S1 of the manuscript and contact the data owners for
the data.

Code required to analyze the diversity data has been uploaded, but note
that the results wont be the same with the manuscript as some private
datasets were not released.

Three data files are available: cleaned_alpha.csv, cleaned_beta.csv, and
cleaned_genus.csv

## Code

In the code folder, there are 7 scripts. Four are crucial for the actual analyses. Please run the scripts in
numerical order

| Script            | Description                                                                 |
|-------------------|-----------------------------------------------------------------------------|
| 1_Alpha           | Alpha diversity analyses                                                    |
| 2_Beta_and_Gamma  | Beta and gamma diversity analyses                                           |
| 3_Diversity_graph | Making graphs presenting diversity changes across land uses (Fig.1 and 2)   |
| 4_Genus           | Genus analyses (including graphs)                                           |
| cleaning_before_release | Removing private datasets and clean data. This doesn't need to be run before running script 1-4 |
| sim               | Simulations about mixed models (Section S5)                                 |
| get_R2            | get R2 from spaMM . Based on the r2() function in the R package performance |

## Alpha diversity analyses

Please use the Cleaned_alpha.csv for these analyses.

| Variable               | Description                                                |
|------------------------|------------------------------------------------------------|
| studyID                | Study identity                                             |
| Lat                    | Latitude                                                   |
| Long                   | Longitude                                                  |
| Exclude.AM.Est.Shannon | Estimated Shannon diversity after excluding Apis mellifera |
| qPD.S                  | Estimated phylogenetic diversity when q (Hill number) = 1  |
| Agr1                   | Agricultural site = 1, non-agricultural site = 0           |
| LC.authoro.corase      | Land use classified based on authors’ descriptions         |
| Nat1                   | Natural site = 1, non-natural site = 0                     |
| Agr1                   | Agricultural site = 1, non-agricultural site = 0           |
| Urb1                   | Urban site = 1, non-urban site = 0                         |
| abundanceMethod        | Sampling methods                                           |
| PC1                    | Climatic PC1 score                                         |
| PC2                    | Climatic PC2 score                                         |
| SC                     | Sampling completeness                                      |

## Beta and gamma diversity analyses

Please use the Cleaned_beta.csv for these analyses. Note that beta
diversity is computed in the script (Gamma/Alpha).

| Variable          | Description                                                |
|-------------------|------------------------------------------------------------|
| study             | Study identity                                             |
| T.gamma           | Taxonomic gamma diversity                                  |
| T.gamma.SC        | Taxonomic gamma diversity sampling completeness            |
| P.gamma           | Phylogenetic gamma diversity                               |
| LC                | Land use (Natural/Agricultural/Urban)                      |
| extent            | Study extent                                               |
| abundanceMethod   | Sampling methods                                           |
| mean_PC1          | Mean climatic PC1 score                                    |
| mean_PC2          | Mean climatic PC2 score                                    |
| SamplingEffort.sd | Standard deviation in sampling efforts within each study   |
| SC.sd             | Sampling completeness standard deviation                   |
| mean_SC           | mean sampling completeness                                 |
| site              | The number of assemblage(or site) in the study             |
| mean_T_alpha      | average rarefied taxonomic alpha diversity in the study    |
| mean_P_alpha      | average rarefied phylogenetic alpha diversity in the study |

## genus analyses

Please use the Cleaned_genus.csv for these analyses.

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
| LC.author.coarse | Land use change (single column)                  |
| orig.value       | same as abundance value                          |
