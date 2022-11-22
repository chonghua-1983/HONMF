# HONMF
HONMF package for manuscript titled "HONMF: integration analysis of multi-omics microbiome data via matrix factorization and hypergraph"

This is an implement of HONMF on multi-omics microbiome data.
Running environment：MATLAB R2019b or later.
The external functions used in the manuscript and JSNMF can be found in the folder external/.

1. HONMF includes the main functions below: <br>

   "OrthNMF.m": HONMF implementation, integrative analysis for multi-omics microbiome data by matrix factorization and hypergraph. <br>
   "OrthNMF_2views.m": an variant of HONMF for addressing to multi-omics data with two modalities. <br>
   "parameter_selection.m": parameters selection for HONMF. <br>
   "parameter_selection_2views.m": parameters selection for OrthNMF_2views.m. <br>
   "sil.py": computing silhouette scores for clustering obtained from each method <br>

2. Reference
