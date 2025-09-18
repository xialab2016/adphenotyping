# Leveraging electronic health records to examine differential clinical outcomes in people with Alzheimer’s Disease

## Identifying AD populations using a novel AD diagnosis phenotyping algorithm - Knowledge driven Online Multimodal Automated Phenotyping (KOMAP)
### See also https://pubmed.ncbi.nlm.nih.gov/37873131/, https://github.com/xinxiong0238/KOMAP, and https://shiny.parse-health.org/KOMAP/

Using an unsupervised phenotyping algorithm, Knowledge driven Online Multimodal Automated Phenotyping (KOMAP), we assigned AD diagnosis status (i.e., probable or possible AD vs not AD) for all patients in the initial AD-EHR data marts. KOMAP is a two step algorithm: (1) KOMAP leverages an online narrative and codified feature search engine (ONCE) powered by multi-source knowledge-graph representation learning to generate a list of informative codified and narrative features relevant to the target clinical concept, i.e., AD. (2) Using the same list of ONCE-selected features, we trained KOMAP independently at UPMC and MGB to account for population heterogeneity. Although KOMAP itself is unsupervised, we evaluated the performance of the AD phenotyping algorithm in predicting AD diagnosis status using gold-standard labels.

## Statistical analysis 

We assessed two clinically relevant indicators of AD decline readily ascertainable from EHR, time to nursing home admission and time to death. We estimated the time-to-outcome using healthcare system-wise covariate-adjusted Cox proportional hazards (PH) models that accounted for the competing risks of nursing home admission and death, stratified by demographic groups. Using patient-level data from both sites, we performed a fixed-effects meta-analysis using inverse variance weighting to estimate the time-to-outcome stratified by demographic groups, adjusting for pooled covariates. 
