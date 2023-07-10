# RElative STrength Of Relative Eigenvectors (RESTORE) 
Descriptions of the measure of RElative STrength Of Relative Eigenvectors (RESTORE) 

## Background
Relative Principal Component Analysis (PCA), namely relative eigenanalysis, deals with the comparison of two covariance matrices. Its mathematical details have been described by Flury.<sup>1,2</sup> The method was introduced into the field of geometric morphometrics by Bookstein and Mitteroecker.<sup>3</sup> An R package named `vcvComp` was developed by Maître and Mitteroecker for implementation of relative PCA with a focus on its applications in geometric morphometrics.<sup>4</sup> 

Relative PCA is implemented through eigendecomposition of $S_B^- S_A$ in `vcvComp`, where $S_A$ and $S_B$ are the covariance matrices of shape variables for two groups, and $S_B^-$ denotes the pseudoinverse of $S_B$. The ratio of generalised variance, defined as the product of all relative eigenvalues, is a scalar summary of the overall relative magnitude of variance between two groups. In the context of geometric morphometrics, a ratio of generalised variance greater than 1 suggests excess shape variance in group A. If groups A and B comprise individuals at juvenile and adult stage of a species, a ratio of generalised variance greater and smaller than 1 would indicate reduced and increased shape variance with development, respectively. 

Large deviations of ratio of generalised variance from 1 suggests variance in one group is large relative to variance in the other group. However, this provides no indication to the impact of relative eigenvectors on observed morphological variations. It is possible to have a large ratio of generalised variance but the actual variance explained by the relative eigenvectors is small. In geometric morphometric analysis of ontogenetic series, this would indicate that although shape variance is strongly developmentally reduced, this would not have a substantial impact on observed morphological variations in the sample. On the other hand, if a large ratio of generalised variance is accompanied by a high proportion of variance explained by the relative eigenvectors, the strongly developmentally constrained shape features would contribute importantly to a converged adult shape (i.e., observed morphological variations would shrink to be primarily lying along the directions of relative eigenvectors).

We propose a measure named RElative STrength Of Relative Eigenvectors (RESTORE) that quantifies the variance of shape along the relative eigenvectors. This measure complements ratio of generalised variance by putting relative eigenvectors into the context of shape variation. Joint interpretation of both measures would “RESTORE” a contextual understanding of the variance between two groups.

## Analytical pipeline for RESTORE score derivation
A schematic illustration of the measure of RESTORE is provided in Fig SX. A step-by-step exposition of the procedures for deriving the RESTORE measure is provided. The R script for calculating the RESTORE measure is provided in Text Sx. We indicate the R script output along with the following descriptions of the mathematical procedures to facilitate understanding of the method and R script.

### Stage 1
Generalized Procrustes Analysis (GPA) is performed based on raw landmark coordinates of both groups (Dimension of $N×pk$, where $N$ is the total sample size in both groups combined, $p$ is the number of landmarks, and $k$ is the dimension of landmarks, which is 2 or 3). This yields GPA-aligned landmark coordinates.

### Stage 2
Principal Component Analysis (PCA) is performed on the GPA-aligned landmark coordinates. This yields a matrix of Principal Component (PC) scores. To ensure numerically stable results, PC scores along the first several eigenvectors were retained for further analysis.<sup>4</sup> We denote the number of dimensions retained by $M$. The eigenvalues ($λ_{1:M}$) represent the magnitude of shape variance along their associated eigenvectors, which are stored in `$EigVal`.

### Stage 3
Relative PCA is performed through eigendecomposition of $S_B^- S_A$, where $S_A$ and $S_B$ is the covariance matrix of the first $M$ PC scores for group A and B, respectively and $S_B^-$ denotes the pseudoinverse of $S_B$. Eigenvectors from eigendecomposition of $S_B^- S_A$ are termed relative eigenvectors. The first and last relative eigenvectors respectively capture shape features with maximum and minimum ratio of variance in group A relative to group B.

Each individual’s PC scores along the first $M$ dimensions are projected onto the directions of relative eigenvectors. This results in an $N×M$ matrix of individuals’ scores along the $M$ relative eigenvectors. Variance of columns of this matrix reflects the magnitude of shape variance along each relative eigenvector $var_{rEigen_{1:M}}$. It should be noted that “shape variance” refers to variation of shape along the first $M$ PCs rather than full shape variance. This is because only the first $M$ PC scores are used in the projection step.

The sum of $var_(rEigen_(1:M))$ divided by shape variance yields the measure of RESTORE, which reflects the proportion of shape variance captured jointly by all $M$ relative eigenvectors. High RESTORE score suggests the relative eigenvectors are generally aligned with the directions of eigenvectors of GPA-aligned landmark coordinates, a circumstance in which the relative eigenvectors would have a stronger impact on observed morphological variations. In contrast, low RESTORE score suggests the relative eigenvectors have a  weak impact on observed morphological variations. We therefore use “Relative strength” to denote the proportion of shape variance captured by relative eigenvectors. A universal benchmark for RESTORE score is irrelevant since high or low RESTORE score should be determined in the context of specific research questions on specific morphological structures.

Two versions of RESTORE measure are calculated, with one having the full shape variance as the denominator (`$RESTORE_All$RESTORE1`) and the other version having shape variance along the first $M$ PCs as the denominator (`$RESTORE_All$RESTORE2`). If `RESTORE1` suggests the proportion of shape variance explained by relative eigenvectors is small, `RESTORE2` can be used to ascertain if this is confounded by limited variation in shape data encoded in $S_A$ and $S_B$ due to the selection of only the first $M$ PCs. Small values of both `RESTORE1` and `RESTORE2` would provide evidence that relative eigenvectors do capture a limited proportion of shape variance, free from bias from using selected shape PCs as input to relative PCA.

### Stage 4
While `RESTORE1` and `RESTORE2` provides two scalar measures of the relative strength of relative eigenvectors, it is also of interest to interrogate each relative eigenvectors with respect to the proportion of total shape variance it captures.

We first subset the $N×M$ matrix of individuals’ scores along the $M$ relative eigenvectors into two submatrices according to whether the associated relative eigenvalues are greater or smaller than 1. The resultant submatrices encode distributions of individuals along the relative eigenvectors featuring shape characteristics with excess variance in group A (`$RESTORE_AExcessVar`) and group B (`$RESTORE_BExcessVar`). Shape variance explained by each relative eigenvector is encoded in `RESTORE_each`. To provide an intuitive understanding of the magnitude of `RESTORE_each`, we also provide the index of eigenvalues ($λ_(1:M)$) closest to each element of `RESTORE_each`, which is stored in `RESTORE_lambdaIdx`. The proportion of shape variance explained by each relative eigenvector relative to full shape variance is encoded in `p_RESTORE_eac`h. 


## Summary
The RESTORE score is essentially assessing whether the directions of relative eigenvectors are aligned with the directions of GPA-aligned landmark coordinates. The score would be high if they are aligned, and the relative eigenvectors would have a strong impact on the observed morphological variations. If the score is low, relative eigenvectors would not strongly impact observed morphological variations, regardless of the magnitude of ratio of generalised variance. RESTORE score is applicable to other research fields where it is necessary to determine alignment of relative eigenvectors with directions of variations of multivariate or high-dimensional data.


## References
1.	Flury B. Some relations between the comparison of covariance matrices and principal component analysis. Comput Stat Data Anal 1983; 1: 97-109.
2.	Flury BN. Analysis of linear combinations with extreme ratios of variance. Journal of the American Statistical Association 1985; 80(392): 915-22.
3.	Bookstein FL, Mitteroecker P. Comparing covariance matrices by relative eigenanalysis, with applications to organismal biology. Evol Biol 2014; 41: 336-50.
4.	Le Maître A, Mitteroecker P. Multivariate comparison of variance in R. Methods Ecol Evol 2019; 10(9): 1380-92.

