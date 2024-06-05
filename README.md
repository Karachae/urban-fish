# urban-fish
Code accompanying "Urbanization is associated with reduced genetic diversity in marine fish populations"

Preprint: https://doi.org/10.1101/2024.02.20.581210 

Data:
* ReProj_FinalUrbanizationdata.csv: Final dataset.

Code:

1.LoopThru_GenDivMeasures.R
* Calculates genetic diversity metrics (gene diversity, allelic richness, population-specific FST, effective population size).
* Raw genetic data not provided.

2.Drivers_Extraction.R
* Extracts, re-projects urbanization metrics (human population density, cumulative human impacts) and merges with previously calculated genetic diverity metrics.
* Creates final dataset.

3.Global_Models.R
* Bayesian generalized linear mixed effects models.
* 32 models: 4 genetic diversity metrics x 2 candidate drivers x 4 buffer sizes.

4.SpatialAutocorr_Check.R
* Tests for spatial autocorrelation in model residuals.
