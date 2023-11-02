library(igraph)
library(mvtnorm)
library(truncnorm)
library(TruncatedNormal)
library(Matrix)

GCut.z = function(z, graph){
  GCut = list()
  for(label in c(0, 1)) {
    where = which(z == label)
    G = igraph::graph.adjacency(graph[where, where], mode = "undirected")
    G = igraph::set.vertex.attribute(G, "node", value = where)
    GCut[[sprintf('%s', label)]] = G
  }
  return(GCut)
}

SWCut = function(z, rho, GCut) {
  ProbEOn = 1 - exp(-rho)
  gList = list()
  gl = 1

  locvec = c()

  for (label in c(0, 1)) {
    G = GCut[[sprintf('%s', label)]]
    E = igraph::E(G)
    V = igraph::V(G)
    vnode = igraph::get.vertex.attribute(G, "node")
    if (length(vnode) == 0){
      warning("Empty subgraph!","\n")
    }
    myz = z[vnode]
    if (length(E) > 0) {
      Eid = stats::rbinom(length(E), 1, ProbEOn[myz[1]+1])
      newG = igraph::delete.edges(G, E[which(Eid == 0)])
      newG = igraph::clusters(newG)
      gList[gl:(gl + newG$no-1)] = lapply(1:newG$no, function(x) {
        igraph::induced.subgraph(graph = G, vids = V[which(newG$membership == x)])
      })

      gl = gl + newG$no
    } else {
      nnode = length(V)
      gList[gl:(gl + nnode-1)] = lapply(1:nnode, function(x) {
        mygraph = igraph::graph.empty(0, directed = FALSE)
        mygraph = igraph::add.vertices(mygraph,nv = 1)
        mygraph = igraph::set.vertex.attribute(mygraph, "node", value = vnode[x])
        return(mygraph)
      })
      gl = gl + nnode
    }
    locvec = c(locvec,vnode)
  }
  return(gList)
}

loglikelihood = function(r, Lambda, z, eta_0, eta_1, sigma2) {
  eta_j_star = eta_0^(1 - z) * eta_1^z
  tmp = Lambda %*% eta_j_star
  loglik = - (sum(r^2) - 2 * sum(r * tmp) + sum(tmp^2)) / (2 * sigma2)
  return(c(loglik))
}

logfullconditional.z = function(r, Lambda, z, eta_0, eta_1, sigma2, pi) {
  loglik = loglikelihood(r, Lambda, z, eta_0, eta_1, sigma2)
  log_fullconditional = loglik + sum(log(pi[z]))
  return(log_fullconditional)
}

update.z = function(SWCutList, r, Lambda, z, eta_0, eta_1, sigma2, pi) {
  updated_z = z
  for (id in 1:length(SWCutList)) {
    current.graph = SWCutList[[id]]
    node.graph = igraph::get.vertex.attribute(current.graph, "node")
    z.graph = z[node.graph]
    newz.graph = (!z[node.graph])*1

    Lambda.graph = as.matrix(Lambda[, node.graph])
    r.graph = r[(rowSums(Lambda.graph) != 0)]

    oldfullconditional = logfullconditional.z(r.graph, Lambda.graph, z.graph, eta_0, eta_1[node.graph], sigma2, pi)
    newfullconditional = logfullconditional.z(r.graph, Lambda.graph, newz.graph, eta_0, eta_1[node.graph], sigma2, pi)

    ratio = exp(newfullconditional - oldfullconditional)
    if (is.infinite(ratio)) {
      update.prob = 1
    } else {
      update.prob = ratio / (1 + ratio)
    }

    is.update = sample(c(0, 1), size = 1, replace = T, prob = c(1 - update.prob, update.prob))
    if (is.update == 1) {
      updated_z[node.graph] = newz.graph
    }
  }
  return(updated_z)
}

loglikelihood2 = function(r, Lambda, z, eta_0, eta_1, sigma2) {
  n = length(r)
  eta_j_star = eta_0^(1 - z) * eta_1^z
  tmp = Lambda %*% eta_j_star
  loglik = n * log(sigma2) / 2 - (sum(r^2) - 2 * sum(r * tmp) + sum(tmp^2)) / (2 * sigma2)
  return(c(loglik))
}

draw_V = function(a_g, b_g, N_g) {
  G = length(N_g)
  sum_N_g = sum(N_g) - cumsum(N_g)

  as = a_g[-G] + N_g[-G]
  bs = b_g[-G] + sum_N_g[-G]

  V = rbeta(n = G - 1, as, bs)
  V = c(V, 1)

  return(V)
}

draw_p = function(V) {
  G = length(V)
  p_g = rep(0, G)
  cum_prod = cumprod(1 - V)
  for (g in 1:(G - 1)) {
    if (g == 1) {
      p_g[g] = V[g]
    } else {
      p_g[g] = cum_prod[g - 1] * V[g]
    }
  }
  p_g[G] = 1 - sum(p_g)

  return(p_g)
}

draw_K_v1 = function(p_vec, m_1, gamma_1, eta_1) {
  k = length(eta_1)
  G = length(p_vec)
  density_mat = matrix(0, nrow = G, ncol = k)
  for (g in 1:G) {
    density_mat[g, ] = dnorm(eta_1, mean = m_1[g], sd = gamma_1[g])
  }
  sample_prob = t(p_vec * density_mat)
  sample_prob = sample_prob / rowSums(sample_prob)
  if (anyNA(sample_prob)) {browser()}
  K = apply(sample_prob, 1, function(x) {sample(G, 1, prob = x)})

  return(K)
}

draw_K = function(p_vec, m_1, gamma_1, eta_1) {
  k = length(eta_1)
  G = length(p_vec)
  log_density_mat = matrix(0, nrow = G, ncol = k)
  for (g in 1:G) {
    log_density_mat[g, ] = dnorm(eta_1, mean = m_1[g], sd = gamma_1[g], log = T)
  }
  log_sample_prob = log(p_vec) + log_density_mat
  #col_min = apply(log_sample_prob, 2, min)
  #log_sample_prob = t(t(log_sample_prob) - col_min)
  #sample_prob = exp(log_sample_prob)
  #sample_prob = sample_prob / rowSums(sample_prob)
  K = apply(log_sample_prob, 2, function(x) {sample_log_prob(size = 1, log_prob_vec = x)})
  if (anyNA(K)) {browser()}

  return(K)
}

sample_log_prob = function(size, log_prob_vec) {
  G = length(log_prob_vec)
  U = runif(G * size)
  gumbel_rv = -log(-log(U))
  values = matrix(gumbel_rv + log_prob_vec, nrow = size, byrow = T)
  values = apply(values, 1, which.max)
  return(values)
}

draw_eta_1_prior = function(K, m_1, gamma_1) {
  k = length(K)
  eta_1 = rnorm(k, mean = m_1[K], sd = sqrt(gamma_1)[K])
  return(eta_1)
}

#' p2fdr
#'
#' function converting posterior inclusion probabilities to FDR
#'
#' @name p2fdr
#' @param posterior_samples a num_samples times k matrix, the posterior samples of latent class labels
#' @return a vector of k, the FDR associated with each of k elements
#' @export

p2fdr = function(posterior_samples) {
  posterior_prob = colMeans(posterior_samples != 0)
  p = posterior_prob
  sorted_p = sort(p, decreasing = T)
  val = cumsum(1 - sorted_p) / (1:length(p))
  order_p = order(p, decreasing = T)
  val[order_p] = val
  return(val)
}

#' FDR_threshold
#'
#' function for generating a threshold on inclusion probabilities based on the input of FDR alpha
#'
#' @name FDR_threshold
#' @param posterior_samples a num_sampels times k matrix, the posterior samples of latent class labels
#' @param alpha a number between 0 and 1, the FDR cutoff
#' @return the threshold on inclusion probabilities based on the input of FDR alpha
#' @export

FDR_threshold = function(posterior_samples, alpha) {
  posterior_prob = colMeans(posterior_samples != 0)
  p = posterior_prob
  sorted_p = sort(p, decreasing = T)
  val = cumsum(1 - sorted_p) / (1:length(p))
  l = which.max(val[val <= alpha])
  threshold = sorted_p[l]
  return(threshold)
}

#' FDR_threshold_from_prob
#'
#' function for generating a threshold on inclusion probabilities based on the input of FDR alpha
#'
#' @name FDR_threshold_from_prob
#' @param posterior_prob a vector of length k, the posterior inclusion probabilities
#' @param alpha a number between 0 and 1, the FDR cutoff
#' @return the threshold on inclusion probabilities based on the input of FDR alpha
#' @export

FDR_threshold_from_prob = function(posterior_prob, alpha) {
  p = posterior_prob
  sorted_p = sort(p, decreasing = T)
  val = cumsum(1 - sorted_p) / (1:length(p))
  l = which.max(val[val <= alpha])
  threshold = sorted_p[l]
  return(threshold)
}

#' BAUM
#'
#' function for Bayesian Analysis for Untargeted Metabolomics data (BAUM)
#'
#' @name BAUM
#' @param r a vector of length p containing feature-specific test statistics
#' @param Q a p times k matrix of matching uncertainty prior information
#' @param g an igraph object of network
#' @param num_samples a positive integer, number of MCMC samples
#' @param m0  an integer, mean for the null component
#' @param G a positive integer, number of clusters for the approximation of Dirichlet process mixture (DPM) of the alternative component
#' @param mu_1 a vector of length G, the prior means of Gaussian distribution in G clusters for the approximation of DPM, i.e., m_g for g = 1, ..., G
#' @param a_1 a number, shape of the Inverse Gamma prior for gamma_0
#' @param a_2 a number, shape of the Inverse Gamma prior for gamma_g
#' @param a_3 a number, shape of the Inverse Gamma prior for sigma^2
#' @param a_4 a number, shape of the Gamma prior for beta_g, i.e., the scale (rate of Gamma) for gamma_g
#' @param a_5 a number, shape of the Inverse Gamma prior for sigma^2_g, i.e., the variance of m_g
#' @param b_1 a number, rate of the Inverse Gamma prior for gamma_0
#' @param b_3 a number, rate of the Inverse Gamma prior for sigma^2
#' @param b_4 a number, rate of the Gamma prior for beta_g, i.e., the scale (rate of Gamma) for gamma_g
#' @param b_5 a number, rate of the Inverse Gamma prior for sigma^2_g, i.e., the variance of m_g
#' @param a_g a vector of length G, Beta prior parameter for p_g for g = 1, .., G
#' @param b_g a vector of length G, Beta prior parameter for p_g for g = 1, .., G
#' @param pi a vector of length 2, prior label distribution; the two elements should be between 0 and 1, and sum up to 1
#' @param rho a vector of 2 positive numbers, global neighborhood parameters
#' @return para_sample, a list of MCMC samples for each parameter, including m_1, beta, gamma_0, gamma_1, sigma_sq, sigma_sq_K, eta_0, eta_1, Lambda, z, N_g and K.
#' @export
#' @examples
#' # this is a toy example to illustrate how to use the BAUM package
#' # data is available under the same folder
#' # metabolites network, features, and potential matchings are a part of the
#' # covid-19 metabolomics data
#' # code and data for this example can be found under the example folder on the github page https://github.com/guoxuan-ma/BAUM
#'
#' # for details, please look at our paper:
#' # Guoxuan Ma, Jian Kang and Tianwei Yu. (2023+). Bayesian Functional Analysis for Untargeted Metabolomics Data with Matching Uncertainty and Very Small Sample Sizes
#'
#' # install package from github
#' devtools::install_github("guoxuan-ma/BAUM")
#'
#' # loading libraries
#' library(BAUM)
#' library(igraph)
#' library(Matrix)
#'
#' # load data
#' load("toy_example.Rdata")
#'
#' # Suppose there are p features and k metabolites. Data should include:
#' # 1. r: the feature-specific summary statistics, should be a vector of length p and can be obtained by statistical tests of user's choice
#' # 2. feature_meta: a binary matrix of dimension p times k, it specifies whether a feature is potentially matched to a metabolite
#' # 3. g: the metabolite network with k vertices
#'
#' # we construct matrix Q, the prior matching confidence
#' # here we assume feature-metabolite has equal matching confidence, which is for illustration but not necessary
#' Q = feature_meta / rowSums(feature_meta)
#'
#' # run BAUM for 1000 MCMC steps
#' set.seed(2023)
#' para_samples = BAUM(r, Q, g, num_samples = 1000)
#'
#' # first 500 steps for burn-in
#' begin = 501
#' end = 1000
#'
#' # FDR threshold at alpha = 0.2 on metabolites posterior inclusion probability
#' FDR_threshold(para_samples$z[begin:end, ], 0.2)
#'
#' # metabolites fdr
#' fdr = p2fdr(para_samples$z[begin:end, ])
#'
#' # prediction on z: whether or not a metabolite is significant
#' pred_z = (colMeans(para_samples$z[begin:end, ]) >= FDR_threshold(para_samples$z[begin:end, ], 0.2)) + 0
#'
#' # plot network structure
#' coords = layout_(g, nicely())
#' plot(g,
#'      layout = coords,
#'      vertex.color = "palegreen3",
#'      vertex.size = 5,
#'      vertex.label = NA,
#'      vertex.frame.color = "white",
#'      edge.width = 0.5)
#'
#' # plot prediction results
#' plot(g,
#'      layout = coords,
#'      vertex.color = c("palegreen3", "salmon")[pred_z + 1],
#'      vertex.size = 5,
#'      vertex.label = NA,
#'      vertex.frame.color = "white",
#'     edge.width = 0.5)
#'
#' # matching uncertainty estimation
#' p = dim(Q)[1]
#' k = dim(Q)[2]
#' pred_Lambda = matrix(0, nrow = p, ncol = k)
#' pred_r_vec = matrix(0, nrow = end - begin + 1, ncol = p)
#' for (s in begin:end) {
#'   pred_Lambda = pred_Lambda + para_samples$Lambda[[s]] / (end - begin + 1)
#'   pred_r_vec[s-begin+1, ] = as.vector(para_samples$Lambda[[s]] %*% ifelse(para_samples$z[s, ],
#'                                                                           para_samples$eta_1[s, ],
#'                                                                           para_samples$eta_0[s]))
#' }
#' pred_Lambda = as.matrix(pred_Lambda)
#'
BAUM = function(r, # test statistics
                Q, # matching uncertainty prior info
                g, # graph
                num_samples = 5000, # number of MCMC samples
                m0 = 0, # mean for the null component
                G = 50, # number of clusters in the alternative component
                mu_1 = 1:50, # mean of m_ki, i.e., mean of Gaussian in G clusters
                a_1 = 10000, # shape of the Inverse Gamma prior for gamma_0
                a_2 = 1000, # shape of the Inverse Gamma prior for gamma_Ki
                a_3 = 2000, # shape of the Inverse Gamma prior for sigma^2
                a_4 = 1000, # shape of the Gamma prior for beta_Ki, the scale (rate of Gamma) for gamma_Ki
                a_5 = 50, # shape of the Inverse Gamma prior for sigma_Ki, the variance of m_ki
                b_1 = 1, # rate of the Inverse Gamma prior for gamma_0
                b_3 = 1000, # rate of the Inverse Gamma prior for sigma^2
                b_4 = 1, # rate of the Gamma prior for beta_Ki, the scale (rate of Gamma) for gamma_Ki
                b_5 = 50, # rate of the Inverse Gamma prior for sigma_Ki, the variance of m_ki
                a_g = rep(10, 50), b_g = rep(10, 50), # Beta prior for p_g's
                pi = c(0.5, 0.5), # prior label distribution
                rho = c(0.1, 0.1) # neighborhood parameter
                ) {
  n = dim(Q)[1]
  k = dim(Q)[2]
  para_samples = list()
  para_samples[["hyperpara"]] = list(m0 = m0, G = G, mu_1 = mu_1,
                                     a_1 = a_1, a_2 = a_2, a_3 = a_3, a_4 = a_4, a_5 = a_5,
                                     b_1 = b_1, b_3 = b_3, b_4 = b_4, b_5 = b_5,
                                     a_g = a_g, b_g = b_g, pi = pi, rho = rho)
  para_samples[['m_1']] = matrix(0, nrow = num_samples, ncol = G)
  para_samples[['beta']] = matrix(0, nrow = num_samples, ncol = G)
  para_samples[['eta_0']] = rep(0, num_samples)
  para_samples[['eta_1']] = matrix(0, ncol = k, nrow = num_samples)
  para_samples[['sigma_sq']] = rep(0, num_samples)
  para_samples[['sigma_sq_K']] = matrix(0, nrow = num_samples, ncol = G)
  para_samples[['gamma_0']] = rep(0, num_samples)
  para_samples[['gamma_1']] = matrix(0, nrow = num_samples, ncol = G)
  para_samples[['Lambda']] = list()
  for (sid in 1:num_samples) {
    para_samples[['Lambda']][[sid]] = Matrix::Matrix(0, nrow = n, ncol = k)
  }
  para_samples[['z']] = matrix(0, ncol = k, nrow = num_samples)
  loglik_list = c()
  para_samples[['N_g']] = matrix(0, nrow = num_samples, ncol = G)
  para_samples[['K']] = matrix(0, nrow = num_samples, ncol = k)
  current_para = list()

  pb = txtProgressBar(min = 0, max = num_samples, style = 3, width = 100, char = "=")
  for (sid in 1:num_samples) {
    # initialization: sample from prior
    if (sid == 1) {
      current_para[['beta']] = rgamma(n = G, a_4, b_4)
      current_para[['gamma_0']] = 1 / rgamma(n = 1, a_1, b_1)
      current_para[['gamma_1']] = 1 / rgamma(n = G, a_2, current_para[['beta']])
      current_para[['sigma_sq']] = 1 / rgamma(n = 1, a_3, b_3)
      current_para[['sigma_sq_K']] = 1 / rgamma(n = G, a_5, b_5)
      current_para[['m_1']] = sort( rnorm(n = G, mu_1, sqrt(current_para[['sigma_sq_K']])) )
      current_para[['N_g']] = rep(k / G, G)
      current_para[['V']] = draw_V(a_g, b_g, current_para[['N_g']])
      current_para[['p']] = draw_p(V = current_para[['V']])
      current_para[['eta_0']] = rtruncnorm(n = 1, a = 0, mean = m0, sd = sqrt(current_para[['gamma_0']]))
      current_para[['K']] = sample(G, size = k, replace = T)
      current_para[['eta_1']] = draw_eta_1_prior(current_para[['K']], current_para[['m_1']], current_para[['gamma_1']])

      current_para[['Lambda']] = matrix(0, nrow = n, ncol = k)
      for (i in 1:n) {
        where = which(Q[i, ] != 0)
        q_i = Q[i, where]
        j = sample(length(where), size = 1, prob = q_i)
        fill.in = rep(0, length(where))
        fill.in[j] = 1
        current_para[['Lambda']][i, where] = fill.in
      }

      current_para[['z']] = rep(1, k)
      # current_para[['z']] = sample(c(0, 1), size = k, replace = T, prob = pi)
      # current_para[['z']] = z
    } else {
      # some preparations
      not_zeros = sum(current_para[['z']] == 1)
      mask = (current_para[['z']] == 0) * 1
      is_z1 = (current_para[['z']] == 1)
      Lambda_J0 = current_para[['Lambda']][, current_para[['z']] == 0]
      Lambda_J1 = current_para[['Lambda']][, current_para[['z']] == 1]
      s = current_para[['Lambda']] %*% mask

      if (not_zeros == (k-1)) {
        Lambda_J0 = matrix(Lambda_J0)
      }

      # update m_1
      sum_eta_1 = rep(0, G)
      have_obs = which(current_para[['N_g']] != 0)
      sum_eta_1[have_obs] = aggregate(current_para[['eta_1']], by = list(current_para[['K']]), FUN = sum)[, 'x']
      w_1 = 1 / (1 / current_para[['sigma_sq_K']] + current_para[['N_g']] / current_para[['gamma_1']])
      v_1 = w_1 * (mu_1 / current_para[['sigma_sq_K']] + sum_eta_1 / current_para[['gamma_1']])
      current_para[['m_1']] = rnorm(n = G, mean = v_1, sd = sqrt(w_1))

      # rearrange clusters according to cluster mean
      ordr = order(current_para[['m_1']])
      current_para[['m_1']] = current_para[['m_1']][ordr]
      current_para[['beta']] = current_para[['beta']][ordr]
      current_para[['gamma_1']] = current_para[['gamma_1']][ordr]
      current_para[['sigma_sq_K']] = current_para[['sigma_sq_K']][ordr]
      current_para[['N_g']] = current_para[['N_g']][ordr]
      current_para[['V']] = current_para[['V']][ordr]
      current_para[['p']] = current_para[['p']][ordr]
      label_map = cbind(1:G, ordr)
      current_para[['K']] = label_map[current_para[['K']], 2, drop = F]

      # mcmc update for eta_0
      sp_sq = 1 / (sum(s^2) / current_para[['sigma_sq']] + 1 / current_para[['gamma_0']])

      if (not_zeros == 1) {
        Lambda_J1 = matrix(Lambda_J1)
        mu_p = sp_sq * (m0 / current_para[['gamma_0']] +
                          sum(s * (r - Lambda_J1 * current_para[['eta_1']][current_para[['z']] == 1])) /  current_para[['sigma_sq']])
      } else {
        mu_p = sp_sq * (m0 / current_para[['gamma_0']] +
                          sum(s * (r - Lambda_J1 %*% current_para[['eta_1']][current_para[['z']] == 1])) /  current_para[['sigma_sq']])
      }
      current_para[['eta_0']] = rtruncnorm(n = 1, a = 0, mean = mu_p, sd = sqrt(sp_sq))

      # mcmc update for eta_1
      if (not_zeros == 0) {
        current_para[['eta_1']] = draw_eta_1_prior(current_para[['K']], current_para[['m_1']], current_para[['gamma_1']])
      } else {
        for (group in 1:G) {
          in_group = (is_z1) & (current_para[['K']] == group)
          num_g = sum(in_group)
          if (num_g == 0) {next}
          Lambda_J1_g = as.matrix(current_para[['Lambda']][, in_group])
          Sigma_p = solve(crossprod(Lambda_J1_g) / current_para[['sigma_sq']] + current_para[['gamma_1']][group] * diag(num_g))
          m_p = tcrossprod(Sigma_p, Lambda_J1_g) %*% (r - current_para[['eta_0']] * rowSums(Lambda_J0)) / current_para[['sigma_sq']] +
            Sigma_p %*% rep(current_para[['m_1']][group], num_g) / current_para[['gamma_1']][group]
          current_para[['eta_1']][in_group] = rmvnorm(n = 1, mean = m_p, sigma = Sigma_p)

        }
        current_para[['eta_1']][mask == 1] = draw_eta_1_prior(current_para[['K']][mask == 1], current_para[['m_1']], current_para[['gamma_1']])
      }

      eta_star = ifelse(current_para[['z']] == 1, current_para[['eta_1']], current_para[['eta_0']])

      # mcmc update for beta
      current_para[['beta']] = rgamma(n = G,
                                      a_2 + a_4,
                                      b_4 + 1 / current_para[['gamma_1']])

      # mcmc update for sigma^2
      arg1 = a_3 + n / 2
      arg2 = b_3 + sum( (r - current_para[['Lambda']] %*% eta_star)^2 ) / 2
      current_para[['sigma_sq']] = 1 / rgamma(n = 1, arg1, arg2)

      # mcmc update for gamma_0
      arg1 = a_1 + 0.5
      arg2 = b_1 + (current_para[['eta_0']] - m0)^2 / 2
      current_para[['gamma_0']] = 1 / rgamma(n = 1, arg1, arg2)

      # mcmc update for gamma_1
      arg1 = a_2 + current_para[['N_g']] / 2
      ss = rep(0, G)
      for (group in 1:G) {
        in_group = (current_para[['K']] == group)
        num_g = sum(in_group)
        if (num_g == 0) {next}
        ss[group] = sum( (current_para[['eta_1']][in_group] - current_para[['m_1']][group])^2 )
      }
      arg2 = current_para[['beta']] + ss / 2
      current_para[['gamma_1']] = 1 / rgamma(n = G, arg1, arg2)

      # update sigma_sq_K
      current_para[['sigma_sq_K']] = 1 / rgamma(n = G,
                                                a_5 + 0.5,
                                                b_5 + (current_para[['m_1']] - mu_1)^2 / 2)

      # update V, p, K and N_g
      current_para[['V']] = draw_V(a_g, b_g, current_para[['N_g']])
      current_para[['p']] = draw_p(V = current_para[['V']])
      current_para[['K']] = draw_K(current_para[['p']], current_para[['m_1']], current_para[['gamma_1']], current_para[['eta_1']])
      current_para[['N_g']] = table(factor(current_para[['K']], levels = 1:G))

      # mcmc update for Lambda
      for (i in 1:n) {
        where = which(Q[i, ] != 0)
        q_i = Q[i, where]
        ll = - (r[i] - eta_star[where])^2 / (2 * current_para[['sigma_sq']])
        log.pos = log(q_i) + ll

        log.ratio = outer(log.pos, log.pos, '-')
        ratio = exp(log.ratio)
        fill.in = rep(0, length(where))
        if (sum(is.infinite(ratio)) > 0) {
          j = which(log.ratio == max(log.ratio), arr.ind = T)[1]
          fill.in[j] = 1
        } else {
          pos = 1 / colSums(ratio)
          j = sample(length(where), size = 1, prob = pos)
          #fill.in = pos
          fill.in[j] = 1
        }

        current_para[['Lambda']][i, where] = fill.in
      }

      # mcmc update for z
      cut = GCut.z(current_para[['z']], g)
      SW = SWCut(current_para[['z']], rho, cut)
      current_para[['z']] = update.z(SW,
                                     r,
                                     current_para[['Lambda']],
                                     current_para[['z']],
                                     current_para[['eta_0']],
                                     current_para[['eta_1']],
                                     current_para[['sigma_sq']],
                                     pi)
    }
    # save samples
    para_samples[['m_1']][sid, ] = current_para[['m_1']]
    para_samples[['beta']][sid, ] = current_para[['beta']]
    para_samples[['gamma_0']][sid] = current_para[['gamma_0']]
    para_samples[['gamma_1']][sid, ] = current_para[['gamma_1']]
    para_samples[['sigma_sq']][sid] = current_para[['sigma_sq']]
    para_samples[['sigma_sq_K']][sid, ] = current_para[['sigma_sq_K']]
    para_samples[['eta_0']][sid] = current_para[['eta_0']]
    para_samples[['eta_1']][sid, ] = current_para[['eta_1']]
    para_samples[['Lambda']][[sid]] = Matrix::Matrix(current_para[['Lambda']])
    para_samples[['z']][sid, ] = current_para[['z']]
    para_samples[['N_g']][sid, ] = current_para[['N_g']]
    para_samples[['K']][sid, ] = current_para[['K']]

    loglik = loglikelihood2(r,
                            current_para[['Lambda']],
                            current_para[['z']],
                            current_para[['eta_0']],
                            current_para[['eta_1']],
                            current_para[['sigma_sq']])
    loglik_list = c(loglik_list, loglik)
    setTxtProgressBar(pb, sid)
  }
  para_samples[['loglik']] = loglik_list
  close(pb)

  return(para_samples)
}

#' generate_hyperpara_grid
#'
#' function to generate a hyperpara grid for sensitivity analysis
#'
#' @name generate_hyperpara_grid
#' @param ratio a vector of 2 positive numbers (r1, r2). Each input parameter is multiplied by r1 and r2 as the value on the parameter grid.
#' @param a_1 a number, shape of the Inverse Gamma prior for gamma_0
#' @param a_2 a number, shape of the Inverse Gamma prior for gamma_g
#' @param a_3 a number, shape of the Inverse Gamma prior for sigma^2
#' @param a_4 a number, shape of the Gamma prior for beta_g, i.e., the scale (rate of Gamma) for gamma_g
#' @param a_5 a number, shape of the Inverse Gamma prior for sigma^2_g, i.e., the variance of m_g
#' @param b_1 a number, rate of the Inverse Gamma prior for gamma_0
#' @param b_3 a number, rate of the Inverse Gamma prior for sigma^2
#' @param b_4 a number, rate of the Gamma prior for beta_g, i.e., the scale (rate of Gamma) for gamma_g
#' @param b_5 a number, rate of the Inverse Gamma prior for sigma^2_g, i.e., the variance of m_g
#' @param beta_pg a number, Beta prior parameter for p_g for g = 1, .., G
#' @return grid, a matrix of parameter grid, each row is a parameter combination
#' @export

generate_hyperpara_grid = function(ratio = c(0.5, 2),
                                   a_1 = 10000, a_2 = 1000, a_3 = 2000, a_4 = 1000, a_5 = 50,
                                   b_1 = 1, b_3 = 1000, b_4 = 1, b_5 = 50, beta_pg = 10) {
  a_1_list = a_1 * ratio
  a_2_list = a_2 * ratio
  a_3_list = a_3 * ratio
  a_4_list = a_4 * ratio
  a_5_list = a_5 * ratio
  b_1_list = b_1
  b_3_list = b_3 * ratio
  b_4_list = b_4
  b_5_list = b_5 * ratio
  beta_pg_list = beta_pg * ratio
  grid = expand.grid(a_1_list, a_2_list, a_3_list, a_4_list, a_5_list, b_1_list, b_3_list, b_4_list, b_5_list, beta_pg_list)
  colnames(grid) = c("a_1", "a_2", "a_3", "a_4", "a_5", "b_1", "b_3", "b_4", "b_5", "beta_pg")
  return(grid)
}
