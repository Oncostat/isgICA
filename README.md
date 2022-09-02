# isgICA
Infinite sparse graphical independent component analysis (Variational Inference)

These codes are the online material for the article:

Rincourt SL, Michiels S, Drubay D. Complex Disease Individual Molecular Characterization Using Infinite Sparse Graphical Independent Component Analysis. Cancer Inform. 2022 Jul 15;21:11769351221105776. doi: [10.1177/11769351221105776](https://doi.org/10.1177%2F11769351221105776). PMID: 35860346; PMCID: PMC9290103.

This article and this code are distributed under the terms of the Creative Commons Attribution-NonCommercial 4.0 License, please cite this article if you used this code.

Correspondance author : `damien.drubay@gustaveroussy.fr`

## Purpose 

Identifying individual mechanisms involved in complex diseases, such as cancer, is essential for precision medicine. Their characterization is particularly challenging due to the unknown relationships of high-dimensional omics data and their inter-patient heterogeneity. We propose to model individual gene expression as a combination of unobserved molecular mechanisms (molecular components) that may differ between the individuals. Considering a baseline molecular profile common to all individuals, these molecular components may represent molecular pathways differing from the population background. 

We defined an infinite sparse graphical independent component analysis (isgICA) to identify these independent molecular components. This model relies on double sparsity: the source matrix sparseness defines the subset of genes involved in each molecular component, whereas the weight matrix sparseness identifies the subset of molecular components associated with each patient. As the number of molecular components is unknown but likely high, we simultaneously inferred it with the weight matrix sparseness using the beta-Bernoulli process (BBP). 

This proposed algorithm provides an insight into individual molecular heterogeneity to better understand complex disease mechanisms. 
We present two folders: one simulation study and one application in early breast cancer.


## How to use ?

### Prerequisites


The R setup is available on the CRAN website (https://cran.r-project.org/).

The following extra-packages are required:
  
  •	zeallot
•	glue
•	GPfit
•	magrittr
•	ggplot2
•	gridExtra
•	reshape2
•	ParBayesianOptimization
•	doParallel
•	compiler
•	class
•	RcppHungarian
•	psych

### Run  

#### Pre-processing

The isgICA model takes in input normalized and center covariables (in row the covariables and in column the individuals).

#### Application

The code of the isgICA model are in the `Rmodel_isgICA.r` file. The data used is in the Breast folder and the simulation study is in the simulation folder. 
To launch the codes, simply set the path in the launch files.



