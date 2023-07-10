# ##############################################################################
# RElative STrength Of Relative Eigenvectors (RESTORE)
# ##############################################################################
# Estimates shape variance captured by each relative eigenvector and all 
# relative eigenvectors combined.

# Argument
# @lmGpA A matrix of N_1 * (p * k), where N_1 is the number of observations in
# group A. p is the number of landmarks and K is the dimension, which is 2 or 3.
# @lmGpB A matrix of N_2 * (p * k), where N_2 is the number of observations in
# group B. p is the number of landmarks and K is the dimension, which is 2 or 3.
# @relEigenObjAB Results from relative.eigen(S1, S2), where S1 is covariance
# matrix for group A and S2 is covariance matrix for group B. Usually scores
# along the first M PCs are used as shape variables in each group.

# Output
# \item{EigVal} Shape varaince explained by each of the first M PCs.
# \item{RESTORE_All} RESTORE_All suggests it contains output from
# the RESTORE function for results obtained by combining all relative
# eigenvectors. It is a list of two components. The first component RESTORE1 
# gives the ratio of shape variance along the first M PCs 
# captured by all relative eigenvectors to shape variance. This is the 
# proportion of shape variance captured by all relative eigenvectors, when 
# relative PCA is performed using the first M PCs of shape. Since relative PCA
# is performed only using M PCs, the proportion of shape variance along the 
# first M PCs captured by all relative eigenvectors is also provided, which is
# stored in RESTORE2. Small RESTORE1 suggests that relative eigenvectors capture 
# a small proportion of shape variance. But it may be possible that this results 
# from the limited shape variance captured by the first M PCs. If RESTORE2 is 
# also small, this would suggest the small variance captured by relative 
# eigenvectors is not due to limited shape variance as input for the relative 
# PCA, but rather is indeed due to the fact that the directions along which 
# variance varies between group A and B capture only a small proportion of shape 
# variance.
# \item{RESTORE_AExcessVar} RESTORE_AExcessVar suggests it contains output from
# the RESTORE function for results on the analysis of the part of the relative
# eigenvectors with relative eigenvalues greater than 1. This is equivalent to 
# the relative eigenvectors where group A has excess variance relative to group
# B. This can be claimed in that S1 and S2 are entered in the the function
# relative.eigen in the order of relative.eigen(S1, S2). RESTORE_each is a 
# numeric vector with length equal to the number of relative eigenvalues greater
# than 1. Each of its component is the shape variance along the first M PCs 
# captured by each relative eigenvector. RESTORE_lambdaIdx contains the index of 
# 𝜆 value closest to each element of RESTORE_each. p_RESTORE_each is the ratio 
# of shape variance along the first M PCs captured by each relative eigenvector 
# to shape variance. This is the proportion of shape variance captured by each
# relative eigenvector.
# \item{RESTORE_BExcessVar} RESTORE_BExcessVar suggests it contains output from
# the RESTORE function for results on the analysis of the part of the relative
# eigenvectors with relative eigenvalues smaller than 1. This is equivalent to 
# the relative eigenvectors where group B has excess variance relative to group
# A. This can be claimed in that S1 and S2 are entered in the the function
# relative.eigen in the order of relative.eigen(S1, S2). RESTORE_each is a 
# numeric vector with length equal to the number of relative eigenvalues smaller
# than 1. Each of its component is the shape variance along the first M PCs 
# captured by each relative eigenvector. RESTORE_lambdaIdx contains the index of 
# 𝜆 value closest to each element of RESTORE_each. p_RESTORE_each is the ratio 
# of shape variance along the first M PCs captured by each relative eigenvector 
# to shape variance. This is the proportion of shape variance captured by each
# relative eigenvector.

# ##############################################################################
RESTORE <- function(lmGpA, lmGpB, relEigenObjAB) {
  
  # PCA of landmark coordinates combining both groups.
  coords <- as.matrix(rbind(lmGpA, lmGpB))
  GpABpca <- prcomp(coords, center = F, scale = F)    
  # ('center' and 'scale' set to F since the covariance matrics are calculated 
  # per group in the cov.group function in vcvComp without centering or 
  # scaling the data a priori, which are then directly used in relative.eigen
  # for eigen decomposition.)
  # Extract PC scores.
  pcScoreMat <- GpABpca$x
  # Select only scores along the first several eigenvectors, with the number
  # determined by the number of relative eigenvectors.
  pcScoreMatRedu <- pcScoreMat[, 1:ncol(relEigenObjAB$relVectors)]
  # (relEigenObjAB and relEigenObjBA have the same number of relative 
  # eigenvectors.)

  # Get the matrix of relative eigenvectors.
  relVec <- relEigenObjAB$relVectors 

  # Project each individuals' PC scores along the selected eigenvectors onto 
  # relative eigenvectors to obtain their scores along relative eigenvectors.
  relVecScore <- pcScoreMatRedu %*% relVec    

  # Calculate variance of scores along relative eigenvectors. This reflects
  # variance of simplied shape (only several eigectors are selected) along the
  # direction of each relative eigenvector.
  RESTORE_each <- apply(relVecScore, 2, var)    
  # Obtain variance of landmark data along each eigenvector.
  EigVal <- (GpABpca$sdev[1:ncol(relVec)])^2

  # For each element of RESTORE_each, find its closest value in all 
  # eigenvalues, not just EigVal. The length of the output should be equal 
  # to the length of RESTORE_each.
  RESTORE_lambdaIdx <- sapply(RESTORE_each, 
                          function(x) (which.min(abs((x - GpABpca$sdev^2)))))


  # Proportion of total shape variance explained by each relative eigenvector.    
  p_RESTORE_each <- RESTORE_each/sum((GpABpca$sdev)^2)
  RESTORE2_each <- RESTORE_each/sum(GpABpca$sdev[1:ncol(relVec)]^2)

  out <- vector(mode = "list", length = 4)
  out[[1]] = EigVal
  out[[2]]$RESTORE1 = sum(p_RESTORE_each)
  out[[2]]$RESTORE2 = sum(RESTORE2_each)
  out[[3]]$RESTORE_each = RESTORE_each[relEigenObjAB$relValues > 1]
  out[[3]]$RESTORE_lambdaIdx = RESTORE_lambdaIdx[relEigenObjAB$relValues > 1]
  out[[3]]$p_RESTORE_each = p_RESTORE_each[relEigenObjAB$relValues > 1]
  out[[4]]$RESTORE_each = rev(RESTORE_each[relEigenObjAB$relValues < 1])
  out[[4]]$RESTORE_lambdaIdx = rev(RESTORE_lambdaIdx[relEigenObjAB$relValues < 1])
  out[[4]]$p_RESTORE_each = rev(p_RESTORE_each[relEigenObjAB$relValues < 1])


  names(out) <- c("EigVal", "RESTORE_All", "RESTORE_AExcessVar", "RESTORE_BExcessVar")
  return(out)

}


# # Example
# library("vcvComp")
# data("Tropheus")

# outliers <- c(18, 56, 155, 351, 624)
# Tropheus.IK <- Tropheus[- outliers, ]
# # Sample reduced to six populations
# Tropheus.IK <- subset(Tropheus.IK, subset = POP.ID %in% levels(POP.ID)[1:6])
# Tropheus.IK$POP.ID <- factor(Tropheus.IK$POP.ID)
# # New variable combining population and sex
# Tropheus.IK$SexPop <- paste(Tropheus.IK$POP.ID, Tropheus.IK$Sex, sep = "_")
# Tropheus.IK$SexPop <- as.factor(Tropheus.IK$SexPop)

# PHEN <- as.matrix(Tropheus.IK[which(names(Tropheus.IK) == "X1"):
# which(names(Tropheus.IK) == "Y19")])
# rownames(PHEN) <- Tropheus.IK$List_TropheusData_ID

# library("geomorph")
# # conversion matrix -> array (19 landmarks, 2 dimensions)
# PHEN_array <- arrayspecs(PHEN, p = 19, k = 2)
# # Procrustes superimposition
# phen.gpa <- gpagen(PHEN_array, print.progress = FALSE)
# # conversion array -> matrix of Procrustes coordinates
# proc.coord <- two.d.array(phen.gpa$coords)
# colnames(proc.coord) <- colnames(PHEN)

# phen.pca <- prcomp(proc.coord, rank. = 5, tol = sqrt(.Machine$double.eps))
# pc.scores <- phen.pca$x

# S.phen.pooled <- cov.group(pc.scores, groups = Tropheus.IK$POP.ID, sex = Tropheus.IK$Sex)
# relEigen.a1s5 <- relative.eigen(S.phen.pooled[, , "IKA1"], S.phen.pooled[, , "IKS5"])
# # RESTORE
# lmIKA1 <- PHEN[which(Tropheus.IK$POP.ID == "IKA1"),]
# lmIKS5 <- PHEN[which(Tropheus.IK$POP.ID == "IKS5"),]
# strg.a1s5 <- RESTORE(lmIKA1, lmIKS5, relEigen.a1s5)
# ##############################################################################