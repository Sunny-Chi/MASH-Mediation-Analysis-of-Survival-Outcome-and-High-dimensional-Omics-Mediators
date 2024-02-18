# MASH: Mediation Analysis of Survival Outcome and High-dimensional Omics Mediators with Application to Complex Diseases
## Description
Environmental exposures such as cigarette smoking influence health outcomes through intermediate molecular phenotypes, such as the methylome, transcriptome, and metabolome. Mediation analysis is a useful tool for investigating the role of potentially high-dimensional intermediate phenotypes in the relationship between environmental exposures and health outcomes. However, little work has been done on mediation analysis when the mediators are high-dimensional and the outcome is a survival endpoint, and none of it has provided a robust measure of total mediation effect. To this end, we propose an estimation procedure for \underline{M}ediation \underline{A}nalysis of \underline{S}urvival outcome and \underline{H}igh-dimensional omics mediators (MASH) based on a second-moment-based measure of total mediation effect for survival data analogous to the $R^2$ measure in a linear model. In addition, we propose a three-step mediator selection procedure to mitigate potential bias induced by non-mediators. Extensive simulations showed good performance of MASH in estimating the total mediation effect and identifying true mediators. By applying MASH to the metabolomics data of 1919 subjects in the Framingham Heart Study, we identified five metabolites as mediators of the effect of cigarette smoking on coronary heart disease risk (total mediation effect, 51.1\%) and two metabolites as mediators between smoking and risk of cancer (total mediation effect, 50.7\%). Application of MASH to a diffuse large B-cell lymphoma genomics data set identified copy-number variations for eight genes as mediators between the baseline International Prognostic Index score and overall survival.
## Installation

To install this package from GitHub, use the following R command:

```r
# install.packages("devtools") # Uncomment if you haven't installed devtools package yet
devtools::install_github("Sunny-Chi/MASH-Mediation-Analysis-of-Survival-Outcome-and-High-dimensional-Omics-Mediators")

