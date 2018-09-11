`%+%` <- function(a, b) paste0(a, b)
# To install R package from a specific github branch:
# devtools::install_github('osofr/simcausal', ref = "master", build_vignettes = FALSE)
# devtools::install_github('osofr/tmlenet', ref = "master", build_vignettes = FALSE)
options(width=120)
require("igraph")
require("simcausal")
require("tmlenet")
require("matrixStats")
require("assertthat")
require("foreach")
source("helper_plotting.R") # source helper plotting functions

# --------------------------------------------------------------------------------------------
# Main sim parameters:
# --------------------------------------------------------------------------------------------
# setwd("/Users/olegsofrygin/GoogleDrive/Network_TMLE/Betsy_sims/")
# bootstrap.var <- TRUE
bootstrap.var <- FALSE
n.bootstrap <- 2
# nsim <- 10
nsim <- 1000
n_replpsi0 <- 1000
# n_replpsi0 <- 5
# simulation sample sizes:
# nsamp_vec <- 10000
nsamp_vec <- c(500, 1000, 10000)

# simulation network sizes (maximum number of friends for each unit):
# K_vec <- seq(5, 30, by=5)
# K_vec <- c(5)
# trimKmax_vec <- list(NULL, 50, 45, 40, 35, 30, 25, 20, 15, 10)
# trimKmax_vec <- list(NULL, 50, 10)
trimKmax_vec <- list(NULL)

# --------------------------------------------------------------------------------------------
# Parallel foreach with "doRedis" package (Redis database needs to be installed beforehand)
# --------------------------------------------------------------------------------------------
outdir <- "./"
# outdir <- "./Rdata_gym_boot/"
# outdir <- "./Rdata_tmp/"
require("itertools")
# require("doRedis")
# registerDoRedis("jobs", password = "JFEFlfki249fkjsk2~.<+JFEFl;") # the password is specific to my redis set-up
# registerDoRedis("jobs") # the usual way to register paralle/l backend

# --------------------------------------------------------------------------------------------
#### to stop all workers from master node:
# --------------------------------------------------------------------------------------------
#### 1. STOP ALL JOBS (from R):
# removeQueue("jobs")
#### 2. FLUSH ALL INFO FROM CURRENT REDIS DB (from shell):
# redis-cli -a "JFEFlfki249fkjsk2~.<+JFEFl;" flushall

#### to from restart redis server after flush (generally not needed)
# sudo service redis-server restart
# redis-server -a "JFEFlfki249fkjsk2~.<+JFEFl;" restart # doesn't seem to work

# --------------------------------------------------------------------------------------------
# To start a local or remote redis worker in R
# (the "jobs" backend registered above will automatically communicate with these workers):
# --------------------------------------------------------------------------------------------
# library("doRedis");
# library("parallel");
# startLocalWorkers(n=3,queue="jobs") # the usual way to start 3 workers locally
# my specific set-up for remote server workers (allows parallelizing over many server nodes)
# startLocalWorkers(n=detectCores(),queue="jobs", password="JFEFlfki249fkjsk2~.<+JFEFl;")
# startLocalWorkers(n=detectCores(),queue="jobs",host=master, password="JFEFlfki249fkjsk2~.<+JFEFl;")

# --------------------------------------------------------------------------------------------
# old way with "doParallel" (should still work on local machine)
# --------------------------------------------------------------------------------------------
require(doParallel)
# registerDoParallel(cores = 1)
registerDoParallel(cores=detectCores())
# (set.cores <- getDoParWorkers())

# --------------------------------------------------------------------------------------------
# Plot re-scaled histogram densities vs. their theoretical standard normal limit
# --------------------------------------------------------------------------------------------
make_ggplot_hist <- function(model_scen, simRes_byscen) {
  require("ggplot2")
  nsamp <- simRes_byscen[["Sim_params"]][["nsamp"]]
  K <- simRes_byscen[["Sim_params"]][["K"]]
  psi0 <- simRes_byscen[["psi0"]]
  # mean stats over sims:
  sims_MeanRes_byscen <- simRes_byscen[["sims_MeanRes_byscen"]]
  # raw results for all sims:
  sims_RawRes_byscen <- simRes_byscen[["sims_RawRes_byscen"]]
  # scen1_Raw_ests <- sims_RawRes_byscen[["Qh.corr"]][["est"]]
  # scen2_Raw_ests <- sims_RawRes_byscen[["Q.miss"]][["est"]]
  # scen3_Raw_ests <- sims_RawRes_byscen[["h.miss"]][["est"]]
  sigma2_0 <- sims_MeanRes_byscen[[model_scen]][,"TrueVar"]
  # sigma2_N <- sims_MeanRes_byscen[[model_scen]][,"depVar.Est"]
  all_psi_n <- sims_RawRes_byscen[[model_scen]][["est"]]

  est_names <- colnames(all_psi_n)
  all_psi_n_scaled <- (all_psi_n - psi0)
  # all_psi_n_scaled <- sqrt(nsamp/K) * (all_psi_n - psi0)
  for (est_name in est_names) {
    print(est_name)
    all_psi_n_scaled[, est_name] <- all_psi_n_scaled[, est_name] / sqrt(sigma2_0[est_name])
  }

  # The theoretical limiting distribution
  vu <- seq(min(all_psi_n_scaled)-0.1, max(all_psi_n_scaled)+0.1, by = .001)
  Normal_distr  <- data.frame(x=rep(vu, length(est_names)), y=rep(dnorm(x = vu, mean = 0, sd = 1), length(est_names)))
  Normal_distr$estimator  <- rep(est_names, each=length(vu))
  all_psi_n_scaled <- data.frame(all_psi_n_scaled)
  all_psi_n_scaled_melted <- reshape2::melt(all_psi_n_scaled)
  colnames(all_psi_n_scaled_melted) <- c("estimator", "estimand")
  hist_by_est <- ggplot(data.frame(all_psi_n_scaled_melted), aes(x=estimand)) +
                  geom_histogram(aes(y=..density..), binwidth=0.1, colour="black", fill="white") +
                  geom_line(data=Normal_distr, aes(x=x, y=y), colour = "red") +
                  # facet_grid(estimator ~ .) +
                  facet_grid(. ~ estimator) +
                  # ggtitle("rescaled $(\\psi_n-\\psi_0)$ vs. N(0,1),\nscenario: " %+% simRes_byscen[["simname"]] %+% ", " %+% model_scen) +
                  ggtitle("rescaled $(\\psi_n-\\psi_0)$ vs. N(0,1), scenario: " %+% simRes_byscen[["simname"]] %+% ", " %+% model_scen) +
                  # theme(plot.title = element_text(lineheight=.8, size = 7))
                  theme(plot.title = element_text(lineheight=.8, size = 7),
                        axis.title = element_text(lineheight=.8, size = 7),
                        panel.margin = grid::unit(0.1, "lines"),
                        plot.margin = grid::unit(c(0,0.4,+0.5,+0.4),"lines")
                        )
  return(hist_by_est)
  # load("netsims.res.all.nsims"%+%500%+%".RData")
}

# --------------------------------------------------------------------------------------------
# Specificy the observed data distribution with a DAG object
# --------------------------------------------------------------------------------------------
create_OnetDAG <- function(K = 5, nC = 1, network = "prefattach", trimKmax = NULL) {
  `%+%` <- function(a, b) paste0(a, b)
  require("igraph")
  require("simcausal")
  options(simcausal.verbose=FALSE)
  if (nC>1) print("generating DAG with nC: " %+% nC)

  # --------------------------------------------------------------------------------------------
  # preferential attachment (BA) model (power law deg distribution):
  # --------------------------------------------------------------------------------------------
  generate.igraph.prefattach <- function(n, power, zero.appeal, m, trimKmax, ...) {
    g <- sample_pa(n, power = power, zero.appeal = zero.appeal, m = m)
    g <- as.directed(as.undirected(g, mode = "collapse"), mode = "mutual")
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(g)
    if (!missing(trimKmax)) {
      # print("trimKmax: "); print(trimKmax)
      NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat, trimKmax = trimKmax)
    } else {
      NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat)
    }
    return(NetInd_out$NetInd_k)
  }
  # --------------------------------------------------------------------------------------------
  # small world (Watts-Strogatz network) model:
  # --------------------------------------------------------------------------------------------
  generate.igraph.smallwld <- function(n, dim, nei, p, ...) {
    # if (!is.null(seed)) set.seed(seed)
    g <- sample_smallworld(dim = 1, size = n, nei = nei, p = p, loops = FALSE, multiple = FALSE)
    g <- as.directed(g, mode = c("mutual"))
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(g) # From igraph object to sparse adj. matrix:
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat) # From igraph object to simcausal/tmlenet input (NetInd_k, nF, Kmax):
    # if (!is.null(seed)) set.seed(NULL)
    return(NetInd_out$NetInd_k)
  }
  # --------------------------------------------------------------------------------------------
  # Network sampler(s) from igraph (regular graph model)
  # Generate regular random graphs with same degree for each node
  # K - degree of each node
  # --------------------------------------------------------------------------------------------
  generate.regular <- function(n, K, ...) {
    if (n <= 200) K <- 5
    igraph.reg <- igraph::sample_k_regular(no.of.nodes = n, k = K, directed = TRUE, multiple = FALSE)
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(igraph.reg)
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat)
    return(NetInd_out$NetInd_k)
  }
  # --------------------------------------------------------------------------------------------
  # sample a network of nC independent communities, same network within each community, total sample size is n
  # --------------------------------------------------------------------------------------------
  gen_network_k_reg_indClusters <- function(n, nC = 1, K, ...) {
    if (n <= 200) {
      K <- 5; nC <- 1
    }
    igraph.reg <- igraph::sample_k_regular(no.of.nodes = as.integer(n/nC), k = K, directed = TRUE, multiple = FALSE)
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(igraph.reg)
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat)
    NetInd_mat_fixed <- NetInd_out$NetInd_k
    NetInd_mat_list <- lapply(seq(nC), function(k_comm) NetInd_mat_fixed + (k_comm - 1)*as.integer(n/nC))
    NetInd_mat <- do.call('rbind', NetInd_mat_list)
    return(NetInd_mat)
  }
  # --------------------------------------------------------------------------------------------
  # Defining the DAG object
  # --------------------------------------------------------------------------------------------
  D <- DAG.empty()
  if (network %in% "prefattach") {
    print("simulating pref attach.")
    # NETWORK MODEL USED IN SCENARIO Rdata_gym_noboot_depW_allScen:
    if (!is.null(trimKmax)) {
      D <- D + network("Net.prefattach", netfun = "generate.igraph.prefattach", power = 0.5, zero.appeal = 5, m = 5, trimKmax = trimKmax)
    } else {
      D <- D + network("Net.prefattach", netfun = "generate.igraph.prefattach", power = 0.5, zero.appeal = 5, m = 5)
    }
  } else if (network %in% "smallworld") {
    print("simulating small world")
    D <- D + network("Net", netfun = "generate.igraph.smallwld", dim = 1, nei = 9, p = 0.1) # small world
  } else if (network %in% "regular") {
    print("simulating regular network graph")
    D <- D + network("Net", netfun = "generate.regular", K = K) # regular (lattice) graph
  } else {
    stop("network undefined")
  }
  D <- D +
      node("latWcat", distr = "rcat.b0", probs = c(0.0494, 0.1823, 0.2806, 0.2680,0.1651, 0.0546)) +
      node("latWnorm", distr = "rnorm", mean = 0, sd = 1) +
      node("HUB", distr = "rconst", const = ifelse(nF >= 25, 1, 0)) # is this person a hub?
  D <- D +
      node("W1", distr = "rcat.b1", probs = c(0.0494, 0.1823, 0.2806, 0.2680,0.1651, 0.0546)) +
      node("W2", distr = "rbern", prob = plogis(-0.2)) +
      node("WNoise", distr = "rbern", prob = plogis(-0.4))

  D <- D +
      node("PA", distr = "rbern", prob = W2*0.05 + (1-W2)*0.15) + # Physically active at baseline (depends on W2)
      node("nF.PA", distr = "rconst", const = sum(PA[[1:Kmax]]), replaceNAw0 = TRUE) # number of phys. active friends
  # Define exposure 0/1 as completely random:
  D <- D + node("A", distr = "rbern", prob = 0.25)
  # Defining the network summary measures based on A:
  D <- D + node("sum.net.A", distr = "rconst", const = (sum(A[[1:Kmax]])*(HUB==0) + sum((W1[[1:Kmax]] > 4)*A[[1:Kmax]])*(HUB==1)), replaceNAw0 = TRUE)

  D <- D +
      node("probY", distr = "rconst",
          const = plogis(ifelse(PA == 1,
                  +5 - 15*(nF.PA < 1), # the probability of maintaining gym membership drops if no friends are PA
                  -8.0 + 0.25*A) +
                  +0.5*sum.net.A + 0.25*nF.PA*sum.net.A + 0.5*nF.PA +
                  +0.5*(W1-1) - 0.58*W2 +
                  -0.5*(3.477-1) + 0.58*0.4496 +
                  +0.5*latWnorm + 0.5*sum(latWnorm[[1:Kmax]])),
          replaceNAw0 = TRUE)
  D <- D + node("Y", distr = "rbern", prob = probY)
  D <- set.DAG(D, latent.v = c("latWcat", "latWnorm", "probY"), n.test = 200)
}

# compute_corrs <- function() {
#   D <- create_OnetDAG(K = 18)
#   datO <- sim(D, n = 50000)
#   print("correlations b/ween W2: ")
#   print(cor(datO[,"W2"], datO[,"W2"][Net_mat[,1]]))
#   print("correlations b/ween corrW.Fj: ")
#   print(cor(datO[,"corrW.F1"], datO[,"corrW.F1"][Net_mat[,1]]))
#   print(cor(datO[,"corrW.F2"], datO[,"corrW.F2"][Net_mat[,2]]))
#   print(cor(datO[,"corrW.F3"], datO[,"corrW.F3"][Net_mat[,3]]))
#   print(cor(datO[,"corrW.F4"], datO[,"corrW.F4"][Net_mat[,4]]))
#   print(cor(datO[,"corrW.F5"], datO[,"corrW.F5"][Net_mat[,5]]))
# }

compute_stats_interventions <- function() {
  # require("tmlenet")
  net.seed <- 123345
  D <- create_OnetDAG(K = 20, network="regular")
  # D <- create_OnetDAG(K = 18, network="prefattach", trimKmax = 40)
  # D <- create_OnetDAG(K = 18, network="prefattach")
  # datO <- sim(D, n = 50000, rndseed = net.seed, rndseed.reset.node = "latWcat")
  datO <- sim(D, n = 20000, rndseed = net.seed, rndseed.reset.node = "A")
  head(datO)
  mean(datO$Y)
  #   ID HUB W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 PA nF.PA A sum.net.A Y
  # 1  1   0  4  0      0        1        1        0        1        0  0     0 0         2 0
  # 2  2   0  4  0      1        1        0        0        0        0  0     0 0         0 0
  # 3  3   0  2  0      0        1        1        1        1        0  0     2 1         2 0
  # 4  4   0  3  1      1        1        0        1        0        0  0     1 0         1 0
  # 5  5   0  2  0      0        1        1        1        1        1  0     1 1         1 0
  # 6  6   0  2  1      0        1        0        0        0        0  0     1 0         1 0
  attributes(datO)$netind_cl
  #   NetInd_k: 2552 716 242 147 299 2090 686 207 1078 613 1961 763 1316 ...
  #   nF: 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5  ...

  # net.seed <- NULL
  # D <- create_OnetDAG(K = 18, network="prefattach")
  # D <- create_OnetDAG(K = 18, network="smallworld")

  # simulate data with fixed network:
  # datO <- sim(D, n = 25000, rndseed = net.seed, rndseed.reset.node = "latWcat")
  head(datO)
  head(datO$W1)
  table(datO$W1)/nrow(datO)
  # simulate data with fixed network and baseline covariates
  # datO <- sim(D, n = 50000, rndseed = net.seed, rndseed.reset.node = "A.obs")

  print(mean(datO$A)) # [1] 0.24776
  print(mean(datO$Y)) # [1] 0.1213 (for N=100000)
  print(mean(datO$PA)) # [1] 0.09928
  mean(datO$Y) - mean(datO$PA) # [1] 0.02066
  print(mean(datO$Y[datO$PA==1])) # [1] 0.5926672
  print(mean(datO$Y[datO$PA==0])) # [1] 0.06934452
  Net_mat <- attributes(datO)$netind_cl$NetInd
  nF <- attributes(datO)$netind_cl$nF

  table(nF)/nrow(datO)
  table(datO$W1)/nrow(datO)
  table(datO$nF.PA)/nrow(datO)
  table(datO$sum.net.A)
  table(datO$sum.net.A)/nrow(datO)
  plot(table(nF))
  plot(table(datO$nF.PA))

  # Mean outcome under no intervention:
  aval <- 0
  Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- Dset + action("gstar", nodes = node("A", distr = "rbern", prob = aset), aset = aval)
  datFull <- sim(Dset, actions="gstar", n = 50000, rndseed = 54321)[["gstar"]]
  Net_mat <- attributes(datFull)$netind_cl$NetInd
  nF <- attributes(datFull)$netind_cl$nF
  (psi0_a0 <- mean(datFull$Y)) # [1] 0.06446
  table(datFull$sum.net.A)/nrow(datFull)
  (psi0_a0_2 <- mean(datFull$Y[nF <= 25])) # [1] 0.0634833

  # Mean outcome under no intervention, but keeping the old summary measures of the hubs unchanged (not interventining on the summaries of hubs)
  aval <- 0
  Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- Dset + action("gstar",
    nodes = c(node("A", distr = "rbern", prob = aset),
              node("sum.net.A", distr = "rconst", const = ifelse(HUB==1, sum.net.A.obs, sum(A[[1:Kmax]])), replaceNAw0 = TRUE)),
    aset = aval)
  datFull <- sim(Dset, actions="gstar", n = 50000, rndseed = 54321)[["gstar"]]
  (psi0_HUB_a0 <- mean(datFull$Y)) # [1] 0.07956
  table(datFull$sum.net.A)/nrow(datFull)

  # 10% random coverage
  aval <- 0.1
  Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- Dset + action("gstar", nodes = node("A", distr = "rbern", prob = aset), aset = aval)
  datFull <- sim(Dset, actions="gstar", n = 50000, rndseed = 54321)[["gstar"]]
  Net_mat <- attributes(datFull)$netind_cl$NetInd
  nF <- attributes(datFull)$netind_cl$nF
  (psi0_a0.1 <- mean(datFull$Y)) # [1] 0.08066
  (psi0_a0.1_2 <- mean(datFull$Y[nF <= 25])) # [1] 0.07760298
  print(psi0_a0.1 - psi0_a0) # [1] 0.01414

  # 25% random coverage
  aval <- 0.25
  Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- Dset + action("gstar", nodes = node("A", distr = "rbern", prob = aset), aset = aval)
  datFull <- sim(Dset, actions="gstar", n = 50000, rndseed = 54321)[["gstar"]]
  Net_mat <- attributes(datFull)$netind_cl$NetInd
  nF <- attributes(datFull)$netind_cl$nF
  (psi0_a0.25 <- mean(datFull$Y)) # [1] 0.12096
  (psi0_a0.25_2 <- mean(datFull$Y[nF <= 25])) # [1] 0.1160422
  print(psi0_a0.25 - psi0_a0)  #[1] 0.05444

  # 40% random coverage
  aval <- 0.4
  Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- Dset + action("gstar", nodes = node("A", distr = "rbern", prob = aset), aset = aval)
  datFull <- sim(Dset, actions="gstar", n = 50000, rndseed = 54321)[["gstar"]]
  Net_mat <- attributes(datFull)$netind_cl$NetInd
  nF <- attributes(datFull)$netind_cl$nF
  (psi0_a0.4 <- mean(datFull$Y)) # [1] 0.19468
  table(datFull$sum.net.A)/nrow(datFull)
  (psi0_a0.4_2 <- mean(datFull$Y[nF <= 25])) # [1] 0.168847
  # NEW ATE:
  print(psi0_a0.4 - psi0_a0.1) # 0.09464

  # 40% random coverage, but keep the old summary measures of the hubs unchanged (not interventining on hubs)
  aval <- 0.4
  Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- Dset + action("gstar",
    nodes = c(node("A", distr = "rbern", prob = aset),
              node("sum.net.A", distr = "rconst", const = ifelse(HUB==1, sum.net.A.obs, sum(A[[1:Kmax]])), replaceNAw0 = TRUE)),
    aset = aval)
  datFull <- sim(Dset, actions="gstar", n = 50000, rndseed = 54321)[["gstar"]]
  Net_mat <- attributes(datFull)$netind_cl$NetInd
  nF <- attributes(datFull)$netind_cl$nF
  (psi0_HUB_a0.4 <- mean(datFull$Y)) # [1] 0.19074
  table(datFull$sum.net.A)/nrow(datFull)
  #         0       1       2       3       4       5       6       7       8       9      10      11      12      13      14      15      16      17      18      19
  # 0.02670 0.10358 0.18786 0.19912 0.15732 0.11152 0.07630 0.05184 0.03432 0.02160 0.01414 0.00790 0.00410 0.00184 0.00120 0.00038 0.00010 0.00006 0.00008 0.00004
  print(psi0_HUB_a0.4 - psi0_a0.1) # 0.09464
  (psi0_a0.4_2 <- mean(datFull$Y[nF <= 25])) # [1] 0.168847

  # TREATMENT EFFECT FROM FULL COVERAGE:
  aval <- 1
  Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- Dset + action("gstar", nodes = node("A", distr = "rbern", prob = aset), aset = aval)
  datFull <- sim(Dset, actions="gstar", n = 50000, rndseed = 54321)[["gstar"]]
  Net_mat <- attributes(datFull)$netind_cl$NetInd
  nF <- attributes(datFull)$netind_cl$nF
  (psi0_a1 <- mean(datFull$Y)) # [1] 0.41434
  (psi0_a1_2 <- mean(datFull$Y[nF <= 25])) # [1] 0.4043262
  print(psi0_a1 - psi0_a0) # [1] 0.34782

  # (1) COVER ONLY THE MOST CONNECTED HUBS (with nF > 15 around %13 of the population and nF >=20 about %6)
  # Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- create_OnetDAG(K = 5, nC = 1, network = "smallworld")
  Dset <- Dset + action("gHubs",
    nodes = c(node("A", distr = "rbern", prob = ifelse(nF >= 20, 0.9, ifelse(nF >= 15, 0.40, 0)))))
  datFull <- sim(Dset, actions="gHubs", n = 50000, rndseed = 54321)[["gHubs"]]
  Net_mat <- attributes(datFull)$netind_cl$NetInd
  nF <- attributes(datFull)$netind_cl$nF
  mean(datFull$A) # [1] 0.09698
  (psi0_gHubs <- mean(datFull$Y)) # [1] 0.131
  (psi0_gHubs_2 <- mean(datFull$Y[nF <= 25])) # [1] 0.1139723
  print(psi0_gHubs - psi0_a0) # [1] [1] 0.05644
  print(psi0_gHubs - psi0_a0.1) # [1] 0.0423

  # (1) Same as above but not intervening on the summaries of HUBs
  # Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- create_OnetDAG(K = 5, nC = 1, network = "smallworld")
  Dset <- Dset + action("gHubs_noHUB",
    nodes = c(node("A", distr = "rbern", prob = ifelse(nF >= 20, 0.9, ifelse(nF >= 15, 0.40, 0))),
              node("sum.net.A", distr = "rconst", const = ifelse(HUB==1, sum.net.A.obs, sum(A[[1:Kmax]])), replaceNAw0 = TRUE))
              )
  datFull <- sim(Dset, actions="gHubs_noHUB", n = 50000, rndseed = 54321)[["gHubs_noHUB"]]
  mean(datFull$A) # [1] 0.09698
  (gHubs_noHUB <- mean(datFull$Y)) # [1] 0.12618
  table(datFull$sum.net.A)/nrow(datFull)

  # (2) INCREASE THE NUMBER OF FRIENDS IN PA BY 1
  # Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- create_OnetDAG(K = 18, network="prefattach")
  # Dset <- create_OnetDAG(K = 18, network="smallworld")
  Dset <- Dset +
    action("plus.nF.PA",
      nodes = node("nF.PA", distr = "rconst", const = ifelse(nF<=15, sum(PA[[1:Kmax]])+1, sum(PA[[1:Kmax]])),
                  replaceNAw0 = TRUE)) +
    # for small-world, due to sparsity, also intervene on everyone who has W1==1 or W1==6
    action("plus.nF.PA.alt",
      nodes = node("nF.PA", distr = "rconst", const = ifelse((nF<=15) | (W1<=2), sum(PA[[1:Kmax]])+1, sum(PA[[1:Kmax]])),
                  replaceNAw0 = TRUE))

  # datFull <- sim(Dset, actions="plus.nF.PA", n = 25000,  rndseed = net.seed, rndseed.reset.node = "latWcat")[["plus.nF.PA"]]
  datFull <- sim(Dset, actions="plus.nF.PA.alt", n = 25000,  rndseed = net.seed, rndseed.reset.node = "latWcat")[["plus.nF.PA.alt"]]

  table(datFull[,"nF.PA"])
  table(datO[,"nF.PA"])
  # head(datFull)
  mean(datFull$A) # [1] 0.2492
  psi0_plusnF.PA <- mean(datFull$Y)
  psi0_plusnF.PA # 0.16976
  psi0_plusnF.PA - print(mean(datO$Y))
  print(psi0_plusnF.PA - psi0_a0) # [1] 0.10324
  print(psi0_plusnF.PA - psi0_a0.1) # [1] 0.0891

  f.gstar1.nFPA <- create_f.gstar_nFPA(nFmax = 10, add.nF.PA = 1)
  nF.PA.new <- f.gstar1.nFPA(datFull)
  dat.new <- cbind(datFull, nF.PA.new = nF.PA.new, diff = nF.PA.new - datFull$nF.PA)
  head(dat.new[dat.new[,"diff"]>0, ])
  max(dat.new[dat.new[,"diff"]>0, "nF"])
  min(dat.new[dat.new[,"diff"]>0, "nF"])

  # JOINT INTERVENTION:
  # (1) COVER ONLY THE MOST CONNECTED HUBS with nF > 15
  # (2) INCREASE THE NUMBER OF FRIENDS IN PA BY 1
  Dset <- create_OnetDAG(K = 5, nC = 1)
  Dset <- Dset + action("plus.nF.PA.gHubs",
    nodes = c(node("nF.PA", distr = "rconst", const = ifelse(nF<=20, sum(PA[[1:Kmax]])+1, sum(PA[[1:Kmax]])), replaceNAw0 = TRUE),
              node("A", distr = "rbern", prob = ifelse(nF >= 20, 0.9, ifelse(nF >= 15, 0.40, 0)))
              ))
  datFull <- sim(Dset, actions="plus.nF.PA.gHubs", n = 50000, rndseed = 54321)[["plus.nF.PA.gHubs"]]
  summary(datFull)

  head(datFull)

  head(datFull[40000:50000,])

  mean(datFull$A) # [1] 0.1038
  psi0_plusnF.PA_gHubs <- mean(datFull$Y)
  psi0_plusnF.PA_gHubs # [1] 0.16698
  print(psi0_plusnF.PA_gHubs - psi0_a0) # [1] 0.14702
  print(psi0_plusnF.PA_gHubs - psi0_a0.1) # [1] 0.08632
}

# --------------------------------------------------------------------------------------------
# Evaluate the true value of the parameter
# (defines the action on the DAG object in Det and simulates full data)
# --------------------------------------------------------------------------------------------
eval_psi0 <- function(gstar1.nodes, gstar2.nodes, ATE = TRUE, nfull = 50000, repl = 10, K, nC = 1, network = "prefattach", net.seed = NULL, rndseed.reset.node = NULL, trimKmax = NULL) {
  `%+%` <- function(a, b) paste0(a, b)
  assertthat::assert_that(is.count(repl))
  if (nC>1) print("evaluting psi0 for nC: " %+% nC)
  print("psi0 eval is using seed: " %+% net.seed %+% "; the seed is reset at node: " %+% rndseed.reset.node)

  if (repl == 1) {
    Dset <- create_OnetDAG(K = K, nC = nC, network = network, trimKmax = trimKmax)
    Dset <- Dset + action("gstar1", nodes = gstar1.nodes)
    Dset <- Dset + action("gstar2", nodes = gstar2.nodes)
    datFull <- sim(Dset, actions=c("gstar1", "gstar2"), n = nfull, rndseed = net.seed, rndseed.reset.node = rndseed.reset.node)
    if (ATE) {
      psi0 <- mean(datFull[["gstar1"]]$Y) - mean(datFull[["gstar2"]]$Y)
    } else {
      psi0 <- mean(datFull[["gstar1"]]$Y)
    }
  # parallel evaluation of psi0:
  } else {
    ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
    ex <- c(ex, "gstar1.nodes", "gstar2.nodes", "nfull", "create_OnetDAG");
    psi0_repl <- foreach(t.counter = 1 : repl, .packages = c("simcausal", "igraph"), .combine = "c", .export = ex) %dopar% {
      Dset <- create_OnetDAG(K = K, nC = nC, network = network, trimKmax = trimKmax)
      Dset <- Dset + action("gstar1", nodes = gstar1.nodes)
      Dset <- Dset + action("gstar2", nodes = gstar2.nodes)
      datFull <- sim(Dset, actions=c("gstar1", "gstar2"), n = nfull, rndseed = net.seed, rndseed.reset.node = rndseed.reset.node)
      if (ATE) {
        psi0 <- mean(datFull[["gstar1"]]$Y) - mean(datFull[["gstar2"]]$Y)
      } else {
        psi0 <- mean(datFull[["gstar1"]]$Y)
      }
      return(psi0)
    }
    psi0 <- mean(psi0_repl)
    print("vector of psi0: "); print(psi0_repl)
  }
  print(c("psi0: ", psi0))
  return(psi0)
}

# --------------------------------------------------------------------------------------------
# Create a stochastic intervention function, this is passed to tmlenet
# --------------------------------------------------------------------------------------------
# Returns a function that samples A with probability aval
# create_f.gstar_stoch <- function(aval, ...) {
#   eval(aval)
#   f.gstar <- function(data, ...){
#     A <- data[,"A", with = FALSE]
#     return(rbinom(n = nrow(data), size = 1, prob = aval[1]))
#   }
#   return(f.gstar)
# }
# # nFmin - minimum number of friends that decides if person is treated or not
# create_f.gstar_dynamic <- function(nFmin = 15, probA = 1, ...) {
#   eval(nFmin)
#   f.gstar <- function(data, ...){
#     nF <- data[,"nF", with = FALSE]
#     return(ifelse(nF >= nFmin[1], rbinom(n = nrow(data), size = 1, prob = probA[1]), 0))
#   }
#   return(f.gstar)
# }
# create_f.gstar_dynamic2 <- function(nFmin1 = 20, nF1.probA = 0.9, nFmin2 = 15, nF2.probA = 0.40) {
#   eval(nFmin1); eval(nFmin2); eval(nF1.probA); eval(nF2.probA);
#   f.gstar <- function(data, ...){
#     nF <- data[,"nF", with = FALSE]
#     A.gstar <- vector(mode="integer", length = length(nF))
#     A.gstar[nF >= 20] <- rbinom(n = sum(nF >= 20), size = 1, prob = nF1.probA[1])
#     A.gstar[(nF >= 15) & (nF < 20)] <- rbinom(n = sum((nF >= 15) & (nF < 20)), size = 1, prob = nF2.probA[1])
#     A.gstar[nF < 15] <- 0L
#     return(A.gstar)
#   }
#   return(f.gstar)
# }
# # network intervention on nF.PA
# create_f.gstar_nFPA <- function(nFmax = 10, add.nF.PA = 1) {
#   eval(nFmax); eval(add.nF.PA)
#   f.gstar <- function(data, ...){
#     nF <- data[,"nF", with = FALSE]
#     nF.PA <- data[,"nF.PA", with = FALSE]
#     return(ifelse(nF <= nFmax[1], nF.PA + add.nF.PA[1], nF.PA))
#   }
#   return(f.gstar)
# }


# --------------------------------------------------------------------------------------------
# Run nsim parallel simulations for a specific scenario in Sim_params (K, nsamp, ...)
# Each run constists of 3 model specification scenarios: Q & h correct, Q misspecified, h misspecified
# Results are saved for: tmle, h_iptw (efficient IPTW), gcomp
# --------------------------------------------------------------------------------------------
runsim <- function(nsim, Sim_params, psi0, n_replpsi0) {
  # --------------------------------------------------------------------------------------------
  # one run of tmlenet, output some stats
  # --------------------------------------------------------------------------------------------
  run_tmlenet <- function(datO, NetInd_mat, sW, sA, K, Qform, hforms, ATE = TRUE, intervene1.sA, intervene2.sA, psi0, bin.method, maxNperBin) {
    if (bin.method%in%"equal.mass") {
      tmlenet_options(bin.method = bin.method, maxNperBin = maxNperBin)
    } else {
      tmlenet_options(bin.method = bin.method, nbins = as.integer(nrow(datO)/maxNperBin))
    }
    options(tmlenet.verbose = FALSE)
    # options(tmlenet.verbose = TRUE)

    # datO_input <- datO[,c("W1", "W2", "W3", "A", "Y")]
    datO_input <- datO
    hform.g0 <- hforms[["hform.g0"]]
    hform.gstar <- hforms[["hform.gstar"]]
    if (as.integer(K)!=as.integer(ncol(NetInd_mat))) K <- ncol(NetInd_mat)
    message(paste0("bin.method: ", bin.method))
    if (!ATE) {
      f.gstar.2 <- NULL
      intervene2.sA <- NULL
    }
    res <- tmlenet(data = datO_input, sW = sW, sA = sA,
                  Ynode = "Y",
                  Kmax = K, NETIDmat = NetInd_mat,
                  intervene1.sA = intervene1.sA,
                  intervene2.sA = intervene2.sA,
                  Qform = Qform,
                  hform.g0 = hform.g0,
                  hform.gstar = hform.gstar,
                  optPars = list(
                    bootstrap.var = bootstrap.var,
                    n.bootstrap = n.bootstrap,
                    boot.nodes = boot.nodes,
                    boot.form = boot.form
                    )
                  )

    extract_ests <- function(estimator, res) {
      # psi0 estimates:
      replace_iptw_name <- rownames(res[[estimator]]$estimates) %in% "h_iptw"
      est <- as.vector(res[[estimator]]$estimates)
      names(est) <- rownames(res[[estimator]]$estimates); names(est)[replace_iptw_name] <- "h.iptw" # to not crash latex
      MAD <- abs(est-psi0[[estimator]])

      # IC-based Vars for dependent obs & IC-based CIs:
      depVar.Est <- as.vector(res[[estimator]]$IC.vars); names(depVar.Est) <- rownames(res[[estimator]]$IC.vars); names(depVar.Est)[replace_iptw_name] <- "h.iptw"
      CIs <- res[[estimator]]$IC.CIs; rownames(CIs)[replace_iptw_name] <- "h.iptw"
      CIcov <- (psi0[[estimator]] <= CIs[,2]) & (psi0[[estimator]] >= CIs[,1])
      depCIlen <- as.vector(CIs[,2]-CIs[,1]); names(depCIlen) <- rownames(CIs)
      # CIlow <- as.vector(CIs[,1]); CIhigh <- as.vector(CIs[,2]);
      # names(CIlow) <- names(CIhigh) <- rownames(CIs)

      # IC_based Vars & CIs conditional on W (allow dependent Q)
      condW.Var.Est <- as.vector(res[[estimator]]$condW.IC.vars); names(condW.Var.Est) <- names(depVar.Est)
      condW.CIs <- res[[estimator]]$condW.CIs; rownames(condW.CIs) <- rownames(CIs)
      condW.CIcov <- (psi0[[estimator]] <= condW.CIs[,2]) & (psi0[[estimator]] >= condW.CIs[,1])
      condW.CIlen <- as.vector(condW.CIs[,2]-condW.CIs[,1]); names(condW.CIlen) <- rownames(CIs)

      # IC_based Vars & CIs conditional on W (assume independent Q)
      condW.indepQ.Var.Est <- as.vector(res[[estimator]]$condW.indepQ.IC.vars); names(condW.indepQ.Var.Est) <- names(depVar.Est)
      condW.indepQ.CIs <- res[[estimator]]$condW.indepQ.CIs; rownames(condW.indepQ.CIs) <- rownames(CIs)
      condW.indepQ.CIcov <- (psi0[[estimator]] <= condW.indepQ.CIs[,2]) & (psi0[[estimator]] >= condW.indepQ.CIs[,1])
      condW.indepQ.CIlen <- as.vector(condW.indepQ.CIs[,2]-condW.indepQ.CIs[,1]); names(condW.indepQ.CIlen) <- rownames(CIs)

      # Parametric bootsrap Vars & CIs:
      bootVar.Est <- as.vector(res[[estimator]]$boot.vars); names(bootVar.Est) <- names(depVar.Est)
      bootCIs <- res[[estimator]]$boot.CIs; rownames(bootCIs) <- rownames(CIs)
      bootCIcov <- (psi0[[estimator]] <= bootCIs[,2]) & (psi0[[estimator]] >= bootCIs[,1])
      bootCIlen <- as.vector(bootCIs[,2]-bootCIs[,1]); names(bootCIlen) <- rownames(CIs)

      # iid Vars & CIs:
      iidVar.Est <- as.vector(res[[estimator]]$iid.vars); names(iidVar.Est) <- names(depVar.Est)
      iidCIs <- res[[estimator]]$iid.CIs; rownames(iidCIs) <- rownames(CIs)
      iidCIcov <- (psi0[[estimator]] <= iidCIs[,2]) & (psi0[[estimator]] >= iidCIs[,1])
      iidCIlen <- as.vector(iidCIs[,2]-iidCIs[,1]); names(iidCIlen) <- rownames(iidCIs)

      # if ((rndseed.reset.node %in% "A") || (rndseed.reset.node %in% "A.obs")) {
      #   depVar.Est <- condW.Var.Est
      #   CIcov <- condW.CIcov
      #   depCIlen <- condW.CIlen
      # }
      # auxilary var ests:
      # aux.vars <- res[[estimator]]$aux.vars
      # print("aux.vars:"); print(aux.vars)

      return(list(condW.indepQ.Var.Est = condW.indepQ.Var.Est, condW.Var.Est = condW.Var.Est, depVar.Est = depVar.Est, bootVar.Est = bootVar.Est, iidVar.Est = iidVar.Est,
                  condW.indepQ.CIcov = condW.indepQ.CIcov, condW.CIcov = condW.CIcov, depCIcov = CIcov, bootCIcov = bootCIcov, iidCIcov = iidCIcov,
                  condW.indepQ.CIlen = condW.indepQ.CIlen, condW.CIlen = condW.CIlen, depCIlen = depCIlen, bootCIlen = bootCIlen, iidCIlen = iidCIlen,
                  MAD = MAD, est = est
                  # , aux.vars = aux.vars
                  ))
    }

    res.summary <- lapply(estimator, extract_ests, res)
    names(res.summary) <- estimator
    return(res.summary)
  }

  # --------------------------------------------------------------------------------------------
  # Simulation parameters
  # --------------------------------------------------------------------------------------------
  scen.name <- Sim_params[["name"]]
  nsamp <- Sim_params[["nsamp"]]
  K <- Sim_params[["K"]]
  trimKmax <- Sim_params[["trimKmax"]]

  network <- Sim_params[["network"]]
  net.seed <- Sim_params[["net.seed"]]
  rndseed.reset.node <- Sim_params[["rndseed.reset.node"]]
  bootstrap.var <- Sim_params[["bootstrap.var"]]
  n.bootstrap <- Sim_params[["n.bootstrap"]]
  boot.nodes <- Sim_params[["boot.nodes"]]
  boot.form <- Sim_params[["boot.form"]]

  ATE <- Sim_params[["ATE"]]
  estimator <- Sim_params[["estimator"]]

  intervene1.sA <- Sim_params[["intervene1.sA"]]
  intervene2.sA <- Sim_params[["intervene2.sA"]]
  gstar1.nodes <- Sim_params[["gstar1.nodes"]]
  gstar2.nodes <- Sim_params[["gstar2.nodes"]]

  # types of models to run:
  Qforms <- Sim_params[["Qforms"]]
  hforms <- Sim_params[["hforms"]]
  bin.methods <- Sim_params[["bin.method"]]
  nobsBin_fracs <- Sim_params[["nobsBin_frac"]]
  nComm <- Sim_params[["nComm"]]
  print("running simulation with nComm: "%+%nComm)

  create_Scen_Params <- function(Q.form, Q.form.name) {
      mapply(function(bin.method, nobsBin_frac) {
                    list(name = Q.form.name%+%"_bin." %+% bin.method %+% "_nobsBinf"%+%nobsBin_frac,
                         bin.method = bin.method, nobsBin_frac = nobsBin_frac,
                         Qform = Q.form, hforms =  hforms[["h.corr"]])
                  }, rep(bin.methods, each = length(nobsBin_fracs)), rep(nobsBin_fracs, times = length(bin.methods)), SIMPLIFY = FALSE)
  }
  Scen_params_Qh.corr <- create_Scen_Params(Q.form = Qforms[["Q.corr"]], Q.form.name = "Qh.corr")
  # Scen_params_Q.miss <- create_Scen_Params(Q.form = Qforms[["Q.miss"]], Q.form.name = "Q.miss")
  Scen_params <- c(Scen_params_Qh.corr)
  # Scen_params <- c(Scen_params_Qh.corr, Scen_params_Q.miss)
  scen_names <- unlist(lapply(Scen_params, "[[", "name"))
  names(Scen_params) <- scen_names

  # --------------------------------------------------------------------------------------------
  # Evaluate the true value of the parameter
  # --------------------------------------------------------------------------------------------
  if (missing(psi0)) {
    psi0 <- vector(mode="list", length=length(estimator))
    names(psi0) <- estimator
    for (est.name in estimator) {
      message("evaluating the true parameter value (psi0)...")
      if (!est.name%in%"ATE") {
        psi0_EY.g1 <- eval_psi0(gstar1.nodes = gstar1.nodes, gstar2.nodes = gstar2.nodes,
                      ATE = FALSE, nfull = nsamp, repl = n_replpsi0, K = K, nC = nComm,
                      network = network, net.seed = net.seed, rndseed.reset.node = rndseed.reset.node, trimKmax = trimKmax)
        psi0[[est.name]] <- psi0_EY.g1
      } else {
        psi0_ATE <- eval_psi0(gstar1.nodes = gstar1.nodes, gstar2.nodes = gstar2.nodes,
                        ATE = ATE, nfull = nsamp, repl = n_replpsi0, K = K, nC = nComm,
                        network = network, net.seed = net.seed, rndseed.reset.node = rndseed.reset.node, trimKmax = trimKmax)
        psi0[[est.name]] <- psi0_ATE
      }
    }
    print("psi0: " %+% psi0)
  }

  # --------------------------------------------------------------------------------------------
  # Run a parallel loop that Simulates observed data then runs tmlenet for each scenario in Scen_params
  # --------------------------------------------------------------------------------------------
  ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
  # ex <- c(ex, "rcat.b1.b0"); # print("ex: "); print(ex)
  simRes <- foreach(t.counter = 1 : nsim, .packages = c("tmlenet","simcausal", "igraph"), .export = ex) %dopar% {
    # ------------------------------------------------------------
    options(simcausal.verbose = FALSE)
    # ------------------------------------------------------------
    # define the data generating distribution for the observed data:
    message("calling create_OnetDAG with trimKmax: " %+% trimKmax)
    Dset <- create_OnetDAG(K = K, nC = nComm, network = network, trimKmax = trimKmax)
    datO <- sim(Dset, n = nsamp, rndseed = net.seed, rndseed.reset.node = rndseed.reset.node)
    netind_cl <- attributes(datO)$netind_cl
    NetInd_mat <- attributes(datO)$netind_cl$NetInd

    res_byQ <- vector(mode="list", length=length(Scen_params)); names(res_byQ) <- names(Scen_params)
    for (scen.name in names(Scen_params)) {
      params <- Scen_params[[scen.name]]
      print("scen.name: " %+% scen.name); print("t.counter: "%+%t.counter)
      print(params$Qform); print(params$hforms)

      tmlenet_res <- run_tmlenet(datO = datO, NetInd_mat = NetInd_mat, K = K,
                                  sW = Sim_params[["sW"]], sA = Sim_params[["sA"]],
                                  Qform = params$Qform, hforms = params$hforms,
                                  ATE = ATE,
                                  intervene1.sA = intervene1.sA,
                                  intervene2.sA = intervene2.sA,
                                  psi0 = psi0,
                                  bin.method = params[["bin.method"]],
                                  maxNperBin = as.integer(params[["nobsBin_frac"]]*nsamp))
      # print("tmlenet_res: "); print(tmlenet_res)
      res_byQ[[scen.name]] <- tmlenet_res
    }
    return(res_byQ)
  }

  # --------------------------------------------------------------------------------------------
  # Collect results and take means over sims:
  # --------------------------------------------------------------------------------------------
  # 1. Break up simRes list into lists by the type of parameter (estimator): EY_gstar1, EY_gstar2 or ATE
  # lapply(simRes, function(onesim_simRes) lapply(onesim_simRes))
  estimator <- unique(names(simRes[[1]][[1]]))
  # the following stacked order: sim index -> names(Scen_params) index -> estimator index (EY_gstar1/ATE)
  simRes_byest <- lapply(estimator, function(est.name) lapply(simRes, function(onesimRes) lapply(onesimRes, '[[', est.name)))
  names(simRes_byest) <- estimator

  # 2. Re-order each parameter sublist by results with only one scenario (a list of simN results by scenario turned into a list of scenarios by simN):
  get_Res_byScen <- function(estimator, simRes_byest) {
    simRes <- simRes_byest[[estimator]]

    simRes_byscen <- lapply(names(Scen_params), function(scen) lapply(simRes, '[[', scen)); names(simRes_byscen) <- names(Scen_params)
    # ouput names from run_tmlenet():
    outparams <- names(simRes_byscen[[1]][[1]])
    # raw results as a list of matrices (each matrix column are results for a particular estimator):
    sims_RawRes_byscen <- lapply(simRes_byscen, function(simRes) lapply(outparams, function(outparam) t(sapply(simRes, '[[', outparam))))

    # sims_RawRes_byscen[["Qh.corr_bin.equal.mass_nobsBinf0.005"]][["est"]]
    # sims_RawRes_byscen[["Qh.corr_bin.equal.mass_nobsBinf0.005"]][["est"]][,"h.iptw"]
    # which(sims_RawRes_byscen[["Qh.corr_bin.equal.mass_nobsBinf0.005"]][["est"]][,"h.iptw"]>1)

    # collect mean matrix results over sims in a list (by scenario):
    sims_MeanRes_byscen <- vector(mode="list", length=length(Scen_params)); names(sims_MeanRes_byscen) <- names(Scen_params)
    # AuxVars_byscen <- vector(mode="list", length=length(Scen_params)); names(AuxVars_byscen) <- names(Scen_params)

    for (scen in names(Scen_params)) {
      names(sims_RawRes_byscen[[scen]]) <- outparams
      # by scenario, take by column means for each output parameter returned by run_tmlenet():
      sims_MeanRes <- lapply(sims_RawRes_byscen[[scen]], base::colMeans)
      bias_abs <- abs(sims_MeanRes[["est"]] - psi0[[estimator]]) # absolute bias
      TrueVar <- colVars(sims_RawRes_byscen[[scen]][["est"]]); names(TrueVar) <- names(sims_MeanRes[["est"]]) # empirical Var
      # AuxVars <- sims_MeanRes[["aux.vars"]]
      sims_MeanRes <- sims_MeanRes[!(names(sims_MeanRes) %in% "aux.vars")]
      sims_MeanRes <- c(list(bias_abs=bias_abs, TrueVar=TrueVar), sims_MeanRes)
      sims_MeanRes <- sapply(sims_MeanRes, I)
      sims_MeanRes_byscen[[scen]] <- sims_MeanRes
      # AuxVars <- matrix(AuxVars, ncol = 1)
      # n_knockout <- c(20, 200, 2000)
      # rownames(AuxVars) <- c("factorized_DY_DW", "noID_1_to_" %+% n_knockout)
      # AuxVars_byscen[[scen]] <- AuxVars
    }

    sims_Summary <- list(simname = Sim_params[["name"]],
                        psi0 = psi0[[estimator]], Sim_params = Sim_params,
                        sims_MeanRes_byscen = sims_MeanRes_byscen,
                        sims_RawRes_byscen = sims_RawRes_byscen
                        # ,
                        # AuxVars_byscen = AuxVars_byscen
                        )
    return(sims_Summary)
  }

  sims_Summary_byest <- lapply(estimator, get_Res_byScen, simRes_byest)
  names(sims_Summary_byest) <- estimator
  return(sims_Summary_byest)
}

# ********************************************************************************************
# Main function that sets-up all simulations
# ********************************************************************************************
# TO DO: Need to fix one network across all sims, evaluate psi0 also on the basis of that single network
# TO DO: Obtain the iid asymptotic var est (ignore dependence), this requires modifying tmlenet() to return the unit-specific ICs
# --------------------------------------------------------------------------------------------
run_all_sims <- function() {
  net.seed <- 123345 # net.seed <- NULL
  # rndseed.reset.node <- "latWcat" # to condition on the network (same network is sampled each time)
  # rndseed.reset.node <- "A.obs"  # to condition on the network and ALL baseline covariates (same network & bsl covars sampled each time)
  rndseed.reset.node <- "A"  # to condition on the network and ALL baseline covariates (same network & bsl covars sampled each time)
  nsim <- nsim
  n_replpsi0 <- n_replpsi0
  nsamp_vec <- nsamp_vec # simulation sample sizes:
  # K_vec <- K_vec # simulation network sizes (maximum number of friends for each unit):
  trimKmax_vec <- trimKmax_vec
  bootstrap.var <- bootstrap.var
  n.bootstrap <- n.bootstrap

  # --------------------------------------------------------------------------------------------
  # Common simulation parameters (These parameters are shared across all simulations)
  # --------------------------------------------------------------------------------------------
  # This should give an error (PA is not defined!)
  sW <-  def_sW(W1, W2, WNoise, HUB = ifelse(nF >= 25, 1, 0)) # PA0 = (PA == 0)
  sA <-  def_sA(A, nF.PA = sum(PA[[1:Kmax]]), replaceNAw0 = TRUE) +
         def_sA(A.PAeq0 = A * (PA == 0)) +
         def_sA(nFPAeq0.PAeq1 = (nF.PA < 1) * (PA == 1)) +
         def_sA(sum.net.A = (sum(A[[1:Kmax]])*(HUB==0) + sum((W1[[1:Kmax]] > 4)*A[[1:Kmax]])*(HUB==1)),
                sum.net.A.sum.netPA = sum.net.A*nF.PA,
                replaceNAw0 = TRUE)

# ------------------------------------------------------------------------------------------------------------------------
# INTERVENING DIRECTLY ON SUMMARY MEASURES: (ADD TO TMLENET VIGNETTE)
# ***** TO DO: combine sW, sA args into one input: obs.sW.sA = c(sW, sA)? *****
# ***** (parameteric bootstrap) What happens when def_new_sA includes A and some summary of A (sum.net.A)?
            # Adding arg 'boot.nodes' to resolve this issue: only resamples nodes specified in boot.nodes,
            # then reconstructs all the summaries in sA
# ------------------------------------------------------------------------------------------------------------------------
# (1) All intervention nodes must be named and must match some previously defined summary measure/node name (been previously defined in def_sA).
# (2) The interention nodes/summaries will replace the existing ones.
# (3) The summaries that were part of def_sA and were not re-defined in def_sA.gstar will be  still re-evaluated on the data generated under def_sA.gstar.
# (4) Each intervention nodes/summary can reference the value of its own previously defined node, evaluated under observed data. For example,
# if we had an observed binary exposure data column A (registered with def_sA(A)), the intervention that reverses the value of A from 0 to 1 and
# from 1 to 0 could be simply defined as def_sA.gstar(A = 1 - A);
# (5) All of the observed exposure summaries defined in obs.sW.sA are evaluated in EXACTLY the same order as they were defined.
# Hence all the intervention summaries preserve exactly the same order of evaluation as in obs.sW.sA;
  # Example 1A: Increase the total number of phys-active friends by 1.
  sA_star1a <- def_new_sA(nF.PA = (nF <= 10)*(sum(PA[[1:Kmax]])+1) + (nF > 10)*sum(PA[[1:Kmax]]), replaceNAw0 = TRUE)
  # # Example 1B: Alternative way of defining exactly the same intervention (will over-ride the existing summary nF.PA defined under sA)
  sA_star1b <- def_new_sA(nF.PA = (nF <= 10)*(nF.PA+1) + (nF > 10)*nF.PA)
  # Example 2: Sample A as a stochastic intervention:
  sA_star2 <-  def_new_sA(A = rbinom(n = length(A), size = 1, prob = 0.10))
  # Example 3A: Sample A as a stochastic intervention and don't intervene on the summaries of the HUBS:
  sA_star3a <- def_new_sA(A = rbinom(n = length(A), size = 1, prob = 0.60)) +
               def_new_sA(sum.net.A = ifelse(HUB==1, sum.net.A, sum(A[[1:Kmax]])), replaceNAw0 = TRUE)
  # Example 3B: Equivalent to 3A, but explicititely defining sum.net.A.sum.netPA as well:
  sA_star3b <- def_new_sA(A = rbinom(n = length(A), size = 1, prob = 0.10)) +
               def_new_sA(sum.net.A = ifelse(HUB==1, sum.net.A, sum(A[[1:Kmax]])), replaceNAw0 = TRUE) +
               def_new_sA(sum.net.A.sum.netPA = sum.net.A*nF.PA)

  # ----------------------------------------------------------------------
  # Regression formulas
  # ----------------------------------------------------------------------
  Qforms <- list(Q.corr = "Y ~ nF.PA + A.PAeq0 + nFPAeq0.PAeq1 + sum.net.A + sum.net.A.sum.netPA + PA + W1 + W2")
  # Qforms <- list(Q.corr = "Y ~ nF.PA + A.PAeq0 + nFPAeq0.PAeq1 + sum.net.A + sum.net.A.sum.netPA + PA + W1 + W2",
  #                Q.miss = "Y ~ A + sum.net.A + W1 + W2")
  # This version creates huge outliers for IPTW at N=1000 (the design mat ordering results in -Inf/+Inf coefficients in speed.glm)
  # hform1.g0 <- "A + sum.net.A ~ HUB + nF.PA + nFPAeq0.PAeq1 + PA0"
  # This ordering does the same speed.glm (-> positivity -> HUGE IPTW):
  # hform2.g0 <- "A + sum.net.A ~ HUB + nF.PA + PA0 + nFPAeq0.PAeq1"
  # this ordering works fine (giving biased, but bounded IPTW):
  hform3.g0 <- "A + sum.net.A ~ HUB + PA + nF.PA + nFPAeq0.PAeq1"

  hforms_A <- list(h.corr = list(hform.g0 = hform3.g0,
                                 hform.gstar = hform3.g0))
  # ----------------------------------------------------------------------
  # STOCHASTIC INTERVENTION ON A
  # ----------------------------------------------------------------------
  # causes positivity issues, switching to more reasonable intervention (0.35)
  # new.sA1.stoch.2 <-  def_new_sA(A = rbinom(n = length(A), size = 1, prob = 0.40))
  # gstar1.nodes.stoch.2 <- c(node("A", distr = "rbern", prob = 0.40))
  new.sA1.stoch.2 <-  def_new_sA(A = rbinom(n = length(A), size = 1, prob = 0.35))
  gstar1.nodes.stoch.2 <- c(node("A", distr = "rbern", prob = 0.35))

  new.sA2.stoch <-  def_new_sA(A = rbinom(n = length(A), size = 1, prob = 0.10))
  gstar2.nodes.stoch <- c(node("A", distr = "rbern", prob = 0.10))

  Common_Sim_Params_2 <- list(scen.name = "2pars.A.stoch0.35",
                            ATE = TRUE, estimator = c("EY_gstar1", "ATE"),
                            boot.nodes = c("A"), boot.form = c("PA ~ W1", "A ~ W2"),
                            intervene1.sA = new.sA1.stoch.2, intervene2.sA = new.sA2.stoch,
                            gstar1.nodes = gstar1.nodes.stoch.2, gstar2.nodes = gstar2.nodes.stoch,
                            nComm = 1, K = 10,
                            sW = sW, sA = sA, nobsBin_frac = c(0.005),
                            bin.method = "equal.mass", Qforms = Qforms, hforms = hforms_A)

  # ----------------------------------------------------------------------
  # DYNAMIC/STOCHASTIC INTERVENTION ON A, CONDITIONAL ON THE NUMBER OF FRIENDS (nF)
  # (sets only around 10% of the population to the exposure)
  # ----------------------------------------------------------------------
  # new.sA1.dyn.4 <-  def_new_sA(A = rbinom(n = length(A), size = 1, prob = ifelse(nF >= 20, 0.9, ifelse(nF >= 15, 0.40, 0))))
  # gstar1.nodes.dyn.4 <- c(node("A", distr = "rbern", prob = ifelse(nF >= 20, 0.9, ifelse(nF >= 15, 0.40, 0))))
  new.sA1.dyn.4 <-  def_new_sA(A = rbinom(n = length(A), size = 1, prob = ifelse(nF >= 20, 0.6, ifelse(nF >= 15, 0.40, 0))))
  gstar1.nodes.dyn.4 <- c(node("A", distr = "rbern", prob = ifelse(nF >= 20, 0.6, ifelse(nF >= 15, 0.40, 0))))
  Common_Sim_Params_4 <- list(scen.name = "2pars.A.dynamic",
                            ATE = TRUE, estimator = c("EY_gstar1", "ATE"),
                            boot.nodes = c("A"), boot.form = c("PA ~ W1", "A ~ W2"),
                            intervene1.sA = new.sA1.dyn.4, intervene2.sA = new.sA2.stoch,
                            gstar1.nodes = gstar1.nodes.dyn.4, gstar2.nodes = gstar2.nodes.stoch,
                            nComm = 1, K = 5,
                            sW = sW, sA = sA, nobsBin_frac = c(0.005),
                            bin.method = "equal.mass", Qforms = Qforms, hforms = hforms_A)

  # ----------------------------------------------------------------------
  # INTERVENING ON NETWORK STRUCTURE (nF.PA) - Add a friend who exercizes, while assigning random 10% to A
  # ----------------------------------------------------------------------
  sW.nFPA <-  def_sW(W1, W2, WNoise, HUB = ifelse(nF >= 25, 1, 0)) +
              def_sW(A) +
              def_sW(sum.net.A = sum(A[[1:Kmax]])*(HUB==0) + sum((W1[[1:Kmax]] > 4)*A[[1:Kmax]])*(HUB==1), replaceNAw0 = TRUE)

  sA.nFPA <-  def_sA(nF.PA = sum(PA[[1:Kmax]]), replaceNAw0 = TRUE) +
              def_sA(A.PAeq0 = A * (PA == 0)) +
              def_sA(nFPAeq0.PAeq1 = (nF.PA < 1) * (PA == 1)) +
              def_sA(nF.PA.PAeq0 = nF.PA * (PA == 0)) +
              def_sA(sum.net.A.sum.netPA = sum.net.A * nF.PA * (PA == 0))

  hforms_nFPA <- list(h.corr = list(hform.g0 = "nF.PA ~ HUB + nF", hform.gstar= "nF.PA ~ HUB + nF"))
  # 1A: Increase the total number of phys-active friends by 1.
  # new.sA1.nFPA.6a <- def_new_sA(nF.PA = (nF <= 15)*(sum(PA[[1:Kmax]])+1) + (nF > 15)*sum(PA[[1:Kmax]]), replaceNAw0 = TRUE)
  # 1B: Alternative way of defining exactly the same intervention (will over-ride the existing summary nF.PA defined under sA)
  new.sA1.nFPA.6a <- def_new_sA(nF.PA = (nF <= 15)*(nF.PA + 1) + (nF > 15)*nF.PA)
  gstar1.nodes.nFPA.6a <- c(node("nF.PA", distr = "rconst", const = ifelse(nF <= 15, sum(PA[[1:Kmax]])+1, sum(PA[[1:Kmax]])), replaceNAw0 = TRUE))

  # No intervention on nF.PA:
  new.sA2.nFPA <- def_new_sA(nF.PA = nF.PA)
  gstar2.nodes.nFPA <- c(node("nF.PA", distr = "rconst", const = sum(PA[[1:Kmax]]), replaceNAw0 = TRUE))

  Common_Sim_Params_6_prefattach <- list(scen.name = "2pars.nFPA.6a",
                            ATE = TRUE, estimator = c("EY_gstar1", "ATE"),
                            boot.nodes = c("PA"), # boot.nodes = c("PA","A"),
                            boot.form = c("PA ~ W2"), # boot.form = c("PA ~ W2", "A ~ W2"),
                            intervene1.sA = new.sA1.nFPA.6a, intervene2.sA = new.sA2.nFPA,
                            gstar1.nodes = gstar1.nodes.nFPA.6a, gstar2.nodes = gstar2.nodes.nFPA,
                            nComm = 1, K = 5,
                            sW = sW.nFPA, sA = sA.nFPA,
                            nobsBin_frac = c(0.005),
                            bin.method = "equal.mass", Qforms = Qforms, hforms = hforms_nFPA
                            )

  # for small world network above node will barely change anybody's nF.PA since most observations have nF>15
  # will also intervene on all observation who have W1<=2:
  new.sA1.nFPA.6b <- def_new_sA(nF.PA = ifelse((nF <= 15) | (W1 <= 2), nF.PA + 1, nF.PA))
  gstar1.nodes.nFPA.6b <- c(node("nF.PA", distr = "rconst", const = ifelse((nF <= 15) | (W1 <= 2), sum(PA[[1:Kmax]])+1, sum(PA[[1:Kmax]])), replaceNAw0 = TRUE))
  Common_Sim_Params_6_smwld <- list(scen.name = "2pars.nFPA.6b",
                            ATE = TRUE, estimator = c("EY_gstar1", "ATE"),
                            boot.nodes =c("PA"), # boot.nodes = c("PA","A"),
                            boot.form = c("PA ~ W2"), # boot.form = c("PA ~ W2", "A ~ W2"),
                            intervene1.sA = new.sA1.nFPA.6b, intervene2.sA = new.sA2.nFPA,
                            gstar1.nodes = gstar1.nodes.nFPA.6b, gstar2.nodes = gstar2.nodes.nFPA,
                            nComm = 1, K = 5,
                            sW = sW.nFPA, sA = sA.nFPA,
                            nobsBin_frac = c(0.005),
                            bin.method = "equal.mass", Qforms = Qforms, hforms = hforms_nFPA
                            )

  # ----------------------------------------------------------------------
  # JOINT INTERVENTION ON THE EXPOSURE (A - dynamic) AND THE NETWORK (nF.PA - add 1)
  # New functionality for specifying interventions (add to ?tmlenet):
  # (0) Allowing multivariate interventions, i.e., one can specify Anodes as a vector of names of intervention nodes that can
  # have length more than 1.
  # (1) Allow multivariate interventions (more than one node) defined by either a list of functions
  # (length(f.g.fun)==length(Anodes)), a single function that returns a matrix of dim n by length(Anodes),
  # a matrix input of the same dimension or a vector of length(Anodes) for constant treatment assignments for each Anode
  # (2) Allowing the user to specify different regression formulas for each intervention node in Anodes
  # (3) If only one regression formula is specified, it will be used for both fitting all of the intervention nodes (with a warning)
  # ----------------------------------------------------------------------
  # JOINT regression for nF.PA and A:

  sW.nFPA.A <-  def_sW(W1, W2, WNoise, HUB = ifelse(nF >= 25, 1, 0))

  sA.nFPA.A <-  def_sA(nF.PA = sum(PA[[1:Kmax]]), replaceNAw0 = TRUE) +
                def_sA(nF.PA.PAeq0 = nF.PA * (PA == 0)) +
                def_sA(nFPAeq0.PAeq1 = (nF.PA < 1) * (PA == 1)) +
                def_sA(A, A.PAeq0 = A * (PA == 0)) +
                def_sA(sum.net.A = sum(A[[1:Kmax]])*(HUB==0) + sum((W1[[1:Kmax]] > 4)*A[[1:Kmax]])*(HUB==1), replaceNAw0 = TRUE) +
                def_sA(sum.net.A.sum.netPA = sum.net.A * nF.PA * (PA == 0))

  hforms_nFPA_A <- list(h.corr = list(
    hform.g0 = c(nF.PA = hforms_nFPA[["h.corr"]][["hform.g0"]], A = hforms_A[["h.corr"]][["hform.g0"]]),
    hform.gstar = c(nF.PA = hforms_nFPA[["h.corr"]][["hform.gstar"]], A = hforms_A[["h.corr"]][["hform.gstar"]])))

  # NEW WAY OF DEFINING JOINT INTERVENTIONS ON THE NETWORK (nF.PA - add 1) AND THE EXPOSURE (A - dynamic):
  new.sA1.nFPA_A.8a <- new.sA1.nFPA.6a + new.sA1.dyn.4
  gstar1.nodes.A_nFPA.8a <- c(gstar1.nodes.nFPA.6a, gstar1.nodes.dyn.4)

  # CONTRAST WITH: (1) no intervention on nF.PA & (2) randomly set 10% to exposure
  new.sA2.nFPA_A <- def_new_sA(nF.PA = nF.PA) +
                    def_new_sA(A = rbinom(n = length(A), size = 1, prob = 0.1))
  gstar2.nodes.A_nFPA <- c(gstar2.nodes.nFPA, gstar2.nodes.stoch)

  Common_Sim_Params_8_prefattach <- list(scen.name = "2pars.A_nFPA.8a",
                            ATE = TRUE, estimator = c("EY_gstar1", "ATE"),
                            boot.nodes = c("PA","A"),
                            boot.form = c("PA ~ W2"), # boot.form = c("PA ~ W1", "A ~ W2"), # only fit model for PA, since A is already fit
                            intervene1.sA = new.sA1.nFPA_A.8a, intervene2.sA = new.sA2.nFPA_A,
                            gstar1.nodes = gstar1.nodes.A_nFPA.8a, gstar2.nodes = gstar2.nodes.A_nFPA,
                            nComm = 1, K = 5,
                            sW = sW.nFPA.A, sA = sA.nFPA.A, nobsBin_frac = c(0.005),
                            bin.method = "equal.mass", Qforms = Qforms, hforms = hforms_nFPA_A)

  # for small world network above node will barely change anybody's nF.PA since most observations have nF>15
  # also intervene on all observation who have W1<=2:
  new.sA1.nFPA_A.8b <- new.sA1.nFPA.6b + new.sA1.dyn.4
  gstar1.nodes.A_nFPA.8b <- c(gstar1.nodes.nFPA.6b, gstar1.nodes.dyn.4)
  Common_Sim_Params_8_smwld <- list(scen.name = "2pars.A_nFPA.8b",
                            ATE = TRUE, estimator = c("EY_gstar1", "ATE"),
                            boot.nodes = c("PA","A"),
                            boot.form = c("PA ~ W2"), # boot.form = c("PA ~ W1", "A ~ W2"), # only fit model for PA, since A is already fit
                            intervene1.sA = new.sA1.nFPA_A.8b, intervene2.sA = new.sA2.nFPA_A,
                            gstar1.nodes = gstar1.nodes.A_nFPA.8b, gstar2.nodes = gstar2.nodes.A_nFPA,
                            nComm = 1, K = 5,
                            sW = sW.nFPA.A, sA = sA.nFPA.A, nobsBin_frac = c(0.005),
                            bin.method = "equal.mass", Qforms = Qforms, hforms = hforms_nFPA_A)


  # --------------------------------------------------------------------------------------------
  # Automatically populates parameter list for all combos of nsamp & trimKmax in nsamp_vec & trimKmax_vec
  # --------------------------------------------------------------------------------------------
 create_Sim_Params <- function(Common_Sim_Params, network) {
      mapply(function(nsamp, trimKmax) {
              c(list(name = "Net_" %+% network %+% "_" %+%
                            Common_Sim_Params$scen.name %+%
                            "_N" %+% nsamp %+%
                            "_Kmax_" %+% ifelse(is.null(trimKmax), "none", trimKmax),
                      nsamp = nsamp, trimKmax = trimKmax, network = network,
                      bootstrap.var = bootstrap.var, n.bootstrap = n.bootstrap,
                      net.seed = net.seed, rndseed.reset.node = rndseed.reset.node),
                      Common_Sim_Params)
              }, rep(nsamp_vec, each = length(trimKmax_vec)), rep(trimKmax_vec, times = length(nsamp_vec)), SIMPLIFY = FALSE)
  }

  # create_Sim_Params(Common_Sim_Params_5b, network = network),
 merge_All_Sim_Params <- function(network = "prefattach") {
    if (network %in% "prefattach") {
      Common_Sim_Params_6 <- Common_Sim_Params_6_prefattach
      Common_Sim_Params_8 <- Common_Sim_Params_8_prefattach
    } else if (network %in% "smallworld") {
      Common_Sim_Params_6 <- Common_Sim_Params_6_smwld
      Common_Sim_Params_8 <- Common_Sim_Params_8_smwld
    }
    All_Sims_Params <- c(
                         create_Sim_Params(Common_Sim_Params_2, network = network),
                         create_Sim_Params(Common_Sim_Params_4, network = network),
                         create_Sim_Params(Common_Sim_Params_6, network = network),
                         create_Sim_Params(Common_Sim_Params_8, network = network)
      )
  }

  All_Sims_Params <- c(merge_All_Sim_Params("prefattach"), merge_All_Sim_Params("smallworld"))
  # All_Sims_Params <- c(merge_All_Sim_Params("regular"))
  # All_Sims_Params <- c(merge_All_Sim_Params("smallworld"))

  length(All_Sims_Params)
  str(All_Sims_Params)
  param_names <- lapply(All_Sims_Params, '[[', "name")
  (names(All_Sims_Params) <- param_names)

  simRes_all <- list()
  for (scenario in names(All_Sims_Params)) {
    time_t <- system.time(simRes_byscen <- runsim(nsim = nsim, Sim_params = All_Sims_Params[[scenario]], n_replpsi0 = n_replpsi0))
    print("sim time for one scenario:"); print(time_t)

    estimators <- All_Sims_Params[[scenario]][["estimator"]]
    simRes_byscen[[estimators]][["Sim_params"]][["scen.name"]]

    for (est.name in estimators) {
      ATE_flag <- simRes_byscen[[est.name]][["Sim_params"]]$ATE
      simRes_byscen[[est.name]][["Sim_params"]]$estimator <- est.name

      simRes_byscen[[est.name]][["simname"]] <- est.name %+% "_" %+% simRes_byscen[[est.name]][["simname"]]
      simRes_byscen[[est.name]][["Sim_params"]][["scen.name"]] <- est.name %+% "_" %+% simRes_byscen[[est.name]][["Sim_params"]][["scen.name"]]

      if (!(est.name %in% "ATE") & ATE_flag) {
        simRes_byscen[[est.name]][["Sim_params"]]$ATE <- FALSE
      }
      simRes_all <- c(simRes_all, list(simRes_byscen[[est.name]]))

      name <- simRes_byscen[[est.name]][["simname"]] # simulation scenario name from All_Sims_Params (Nsamp and K) - should match w/ scenario var
      sims_MeanRes_byscen <- simRes_byscen[[est.name]][["sims_MeanRes_byscen"]] # mean result tables by model spec scenario
      # sims_AuxVars_byscen <- simRes_byscen[[est.name]][["AuxVars_byscen"]] # mean result table for various variance estimates
      sims_RawRes_byscen <- simRes_byscen[[est.name]][["sims_RawRes_byscen"]] # raw simulation results by model spec scenario
      model_names <- names(sims_MeanRes_byscen) # names of different model spec scenarios for Q & h:
      print("results for scenario: " %+% name); print(sims_MeanRes_byscen)
      # print("auxilary variance estimates for scenario: " %+% name); print(sims_AuxVars_byscen)
      # --------------------------------------------------------------------------------------------
      # Write sim result table to a text file
      # --------------------------------------------------------------------------------------------
      sink(file=outdir%+%"SimResults.txt", append=T)
      cat("\n#results for scenario: ", name, "\n")
      print(sims_MeanRes_byscen, quote = FALSE)
      # print(sims_AuxVars_byscen, quote = FALSE)
      sink()
      # --------------------------------------------------------------------------------------------
      # create ggplot histograms for each model scenario:
      # --------------------------------------------------------------------------------------------
      ggplothist_byscen <- lapply(model_names, make_ggplot_hist, simRes_byscen = simRes_byscen[[est.name]])
      # hist <- make_ggplot_hist(model_scen = "Qh.corr", simRes_byscen = simRes_byscen)
      names(ggplothist_byscen) <- model_names
      pdf(file=outdir%+%"fig_" %+% scenario %+% "_" %+% est.name %+%".pdf",width=8,height=12) # pdf flip: width=12,height=8
      multiplot(plotlist=ggplothist_byscen, cols=1) # cols=length(ggplothist_byscen))
      dev.off()

      # names(simRes_all) <- names(All_Sims_Params)[1:length(simRes_all)]
      names(simRes_all) <- unlist(lapply(simRes_all, '[[', "simname"))
    }
    # --------------------------------------------------------------------------------------------
    # save intermediate results:
    # --------------------------------------------------------------------------------------------
    save(list=c("simRes_all"), file=outdir%+%"simRes_all.RData")

  }
  save(list=c("simRes_all"), file=outdir%+%"simRes_all.RData")
  # load(file=outdir%+%"simRes_all.RData")
}

run_all_sims()



