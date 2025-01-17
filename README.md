# KOMAP

We used a novel algorithm, Knowledge driven Online Multimodal Automated Phenotyping (KOMAP), to assign AD diagnosis status (probable / possible AD, or not AD) for all patients in the AD data mart. KOMAP is a two-step unsupervised algorithm that enables retrieval of accurate diagnostic labels from EHR data. In the first step, KOMAP leverages an online narrative and codified feature search engine (ONCE) powered by multi-source knowledge graph representation learning to generate a list of informative codified and narrative features relevant to AD. In the second step, we trained the KOMAP algorithm based on the features selected by ONCE. These features include highly predictive surrogate features of AD diagnosis (e.g., main AD PheCode from codified data, main AD CUI from narrative data), and a measure of healthcare utilization (defined as total number of ICD codes and clinical encounters). 

# Disparity in AD-related outcomes

We examine disparity in two clinically meaningful AD-related outcomes, time-to-nursing home admission and time-to-death.
