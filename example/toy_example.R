# this is a toy example to illustrate how to use the BAUM package
# data is available under the same folder
# metabolites network, features, and potential matchings are a part of the
# covid-19 metabolomics data

# for details, please look at our paper:
# Guoxuan Ma, Jian Kang and Tianwei Yu. (2023+). Bayesian Functional Analysis for Untargeted Metabolomics Data with Matching Uncertainty and Very Small Sample Sizes

# install package from github
devtools::install_github("guoxuan-ma/BAUM")

# loading libraries
library(BAUM)
library(igraph)
library(Matrix)

# load data
load("toy_example.Rdata")

# Suppose there are p features and k metabolites. Data should include:
# 1. r: the feature-specific summary statistics, should be a vector of length p and can be obtained by statistical tests of user's choice
# 2. feature_meta: a binary matrix of dimension p times k, it specifies whether a feature is potentially matched to a metabolite
# 3. g: the metabolite network with k vertices

# we construct matrix Q, the prior matching confidence
# here we assume feature-metabolite has equal matching confidence, which is for illustration but not necessary
Q = feature_meta / rowSums(feature_meta)

# run BAUM for 1000 MCMC steps
set.seed(2023)
para_samples = BAUM(r, Q, g, num_samples = 1000)

# first 500 steps for burn-in
begin = 501
end = 1000

# FDR threshold at alpha = 0.2 on metabolites posterior inclusion probability
FDR_threshold(para_samples$z[begin:end, ], 0.2)

# metabolites fdr
fdr = p2fdr(para_samples$z[begin:end, ])

# prediction on z: whether or not a metabolite is significant
pred_z = (colMeans(para_samples$z[begin:end, ]) >= FDR_threshold(para_samples$z[begin:end, ], 0.2)) + 0

# plot network structure
coords = layout_(g, nicely())
plot(g,
     layout = coords,
     vertex.color = "palegreen3",
     vertex.size = 5,
     vertex.label = NA,
     vertex.frame.color = "white",
     edge.width = 0.5)

# plot prediction results
plot(g,
     layout = coords,
     vertex.color = c("palegreen3", "salmon")[pred_z + 1],
     vertex.size = 5,
     vertex.label = NA,
     vertex.frame.color = "white",
     edge.width = 0.5)

# matching uncertainty estimation
p = dim(Q)[1]
k = dim(Q)[2]
pred_Lambda = matrix(0, nrow = p, ncol = k)
pred_r_vec = matrix(0, nrow = end - begin + 1, ncol = p)
for (s in begin:end) {
  pred_Lambda = pred_Lambda + para_samples$Lambda[[s]] / (end - begin + 1)
  pred_r_vec[s-begin+1, ] = as.vector(para_samples$Lambda[[s]] %*% ifelse(para_samples$z[s, ],
                                                                          para_samples$eta_1[s, ],
                                                                          para_samples$eta_0[s]))
}
pred_Lambda = as.matrix(pred_Lambda)
