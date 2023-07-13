# RElative STrength Of Relative Eigenvectors (RESTORE) 
<!-- ${\color{#D41C5C}RE}{lative}$ ${\color{#00B0F0}ST}{rength}$ ${\color{#FF00FF}O}{f}$ ${\color{#9900FF}R}{elative}$ ${\color{#70AD47}E}{igenvectors}$ -->

![RESTOREchart](https://github.com/Patrick-Wen/RESTORE/assets/100295693/b0d31458-6635-4d4a-be93-cf2b896a9cc1)

## An introduction to RElative STrength Of Relative Eigenvectors (RESTORE)

### Background
Relative Principal Component Analysis (PCA), namely relative eigenanalysis, deals with the comparison of two covariance matrices. The method is established by Flury in the 1980s.<sup>1,2</sup> It was introduced into the field of geometric morphometrics by Bookstein and Mitteroecker.<sup>3</sup> An R package named `vcvComp` was developed by Maître and Mitteroecker for implementation of relative PCA with a focus on its applications in geometric morphometrics.<sup>4</sup> 

Relative PCA identifies directions with high contrasts of variance between two groups. The first relative eigenvector finds the direction where group A has maximum excess variance over group B. The last relative eigenvector is the direction where group A has minimum excess variance, which is also interpreted as the direction where group B has maximum excess variance relative to group A. A summary of the relative magnitude of variance between two groups along all relative eigenvectors is given by the ratio of generalised variance. In the context of geometric morphometrics, if group A and B comprise individuals at juvenile and adult stage of a species, a ratio of generalised variance greater and smaller than 1 would indicate reduced and increased shape variance with development, respectively.

It should be noted that directions of high variance contrasts can be different from the directions of variation in the pooled data (i.e., combining data from both groups). If there is little variation of data along the relative eigenvectors, the relative eigenvectors would not effectively explain the pattern of variation in the data, irrespective of the absolute value of ratio of generalised variance. In geometric morphometric analysis of ontogenetic series, this corresponds to a situation where although ratio of generalised variance suggests shape variance is developmentally reduced, the observed morphological variation in the pooled sample is minimally affected. Results from relative PCA should then not be overinterpreted in that relative eigenvectors barely capture any variation in the data.

It is useful to develop a measure that reflects the proportion of data variation explained by relative eigenvectors. We propose RElative STrength Of Relative Eigenvectors (RESTORE) to identify situations where results from relative PCA should be interpreted with caution. RESTORE would “restore” a contextual understanding of the variance between two groups by shifting attention away from a mere focus inside the relative eigenspace to its evaluation in the eigenspace of pooled multivariate data.

### Analytical pipeline for RESTORE score derivation
A schematic illustration of the measure of RESTORE is provided in Fig SX. A step-by-step exposition of the procedures for deriving the RESTORE measure is provided. The R script for calculating the RESTORE score is provided in Text Sx. In the following descriptions of the analytical pipeline, output from the R script is indicated where appropriate to facilitate understanding of the method and R script.

#### Stage 1
Generalized Procrustes Analysis (GPA) is performed based on raw landmark coordinates of both groups (Dimension of $N×pk$, where $N$ is the total sample size in both groups combined, $p$ is the number of landmarks, and $k$ is the dimension of landmarks, which is 2 or 3). This yields GPA-aligned landmark coordinates.

#### Stage 2
Principal Component Analysis (PCA) is performed on the GPA-aligned landmark coordinates. This yields a matrix of Principal Component (PC) scores. To ensure numerically stable results, PC scores along the first several eigenvectors were retained for further analysis.<sup>4</sup> We denote the number of dimensions retained by $M$. The eigenvalues ($λ_{1:M}$) represent the magnitude of shape variance along their respective eigenvectors, which are stored in `$EigVal`.

#### Stage 3
Relative PCA is performed through eigendecomposition of $S_B^- S_A$, where $S_A$ and $S_B$ is the covariance matrix of the first $M$ PC scores for group A and B, respectively and $S_B^-$ denotes the pseudoinverse of $S_B$. Eigenvectors from eigendecomposition of $S_B^- S_A$ are termed relative eigenvectors. The first and last relative eigenvector respectively captures shape features with maximum and minimum ratio of variance in group A relative to group B.

Each individual’s PC scores along the first $M$ dimensions are projected onto the directions of relative eigenvectors. This results in an $N×M$ matrix of individuals’ scores along the $M$ relative eigenvectors. Variance of columns of this matrix reflects the magnitude of shape variance along each relative eigenvector ($var_{rEigen_{1:M}}$). It should be noted that “shape variance” refers to variation of shape along the first $M$ PCs rather than full shape variance. This is because only the first $M$ PC scores are used in the projection step.

The sum of $var_{rEigen_{1:M}}$ divided by shape variance yields the measure of RESTORE, which reflects the proportion of shape variance captured jointly by all $M$ relative eigenvectors. High RESTORE score suggests the relative eigenvectors are generally aligned with the directions of eigenvectors of GPA-aligned landmark coordinates, a circumstance in which the relative eigenvectors would have a relatively strong impact on observed morphological variations. In contrast, low RESTORE score suggests the relative eigenvectors have a weak impact on observed morphological variations. This also explains why this measure is given the name “RElative STrength Of Relative Eigenvectors”.

Two versions of RESTORE measure are calculated, with one having the full shape variance as the denominator (`$RESTORE_All$RESTORE1`) and the other version having shape variance along the first $M$ PCs as the denominator (`$RESTORE_All$RESTORE2`). If `RESTORE1` suggests the proportion of shape variance explained by relative eigenvectors is small, `RESTORE2` can be used to ascertain if this is confounded by limited variation in shape data encoded in $S_A$ and $S_B$ due to the selection of only the first $M$ PCs. Small values of both `RESTORE1` and `RESTORE2` would provide evidence that relative eigenvectors do capture a limited proportion of shape variance, free from bias from using selected shape PCs as input to relative PCA.

#### Stage 4
While `RESTORE1` and `RESTORE2` provides two scalar measures of the relative strength of relative eigenvectors, it is also of interest to interrogate each relative eigenvectors individually to determine the proportion of total shape variance it captures.

We first subset the $N×M$ matrix of individuals’ scores along the $M$ relative eigenvectors into two submatrices according to whether the associated relative eigenvalues are greater or smaller than 1. The resultant submatrices encode distributions of individuals along the relative eigenvectors featuring shape characteristics with excess variance in group A (`$RESTORE_AExcessVar`) and group B (`$RESTORE_BExcessVar`). Shape variance explained by each relative eigenvector is encoded in `RESTORE_each`. To provide an intuitive understanding of the magnitude of `RESTORE_each`, we also provide the index of eigenvalues ($λ_{1:M}$) closest to each element of `RESTORE_each`, which is stored in `RESTORE_lambdaIdx`. The proportion of shape variance explained by each relative eigenvector relative to full shape variance is encoded in `p_RESTORE_each`. 


### Interpretation
The merit of RESTORE score is that a low value would suggest that the difference in covariance structure between two groups is a trivial explanation for the observed variation (i.e., morphological variation in the context of geometric morphometrics) in the pooled data. Results from relative PCA should not be overinterpreted in this condition. A high value of RESTORE score reflects that data are indeed dispersed along the directions of high variance contrasts, suggesting that between-group difference in covariance structure nontrivially explains the observed data variation.

### Summary
The RESTORE score is essentially assessing whether the directions of relative eigenvectors are aligned with the directions of variation of pooled data (GPA-aligned landmark coordinates in the context of geometric morphometrics). Low RESTORE score is of particular interest in that it suggests relative eigenvectors do not strongly impact the observed morphological variations and results from relative PCA should be interpreted with caution.


### References
1.	Flury B. Some relations between the comparison of covariance matrices and principal component analysis. Comput Stat Data Anal 1983; 1: 97-109.
2.	Flury BN. Analysis of linear combinations with extreme ratios of variance. Journal of the American Statistical Association 1985; 80(392): 915-22.
3.	Bookstein FL, Mitteroecker P. Comparing covariance matrices by relative eigenanalysis, with applications to organismal biology. Evol Biol 2014; 41: 336-50.
4.	Le Maître A, Mitteroecker P. Multivariate comparison of variance in R. Methods Ecol Evol 2019; 10(9): 1380-92.

