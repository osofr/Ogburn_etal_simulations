`%+%` <- function(a, b) paste0(a, b)
# To install R package from a specific github branch:
# devtools::install_github('osofr/simcausal', ref = "master", build_vignettes = FALSE)
# devtools::install_github('osofr/tmlenet', ref = "master", build_vignettes = FALSE)
options(width=120)
require("simcausal")
require("tmlenet")
require("matrixStats")
require("assertthat")
source("helper_plotting.R") # source helper plotting functions
# --------------------------------------------------------------------------------------------
# Parallel foreach with "doRedis" package (Redis database needs to be installed beforehand)
# --------------------------------------------------------------------------------------------
# outdir <- "./Rdata/"
# outdir <- "./Rdata_tmp/"
outdir <- "./Rdata_catA_depW_nC1/"
require("itertools")
require("doRedis")
# registerDoRedis("jobs") # the usual way to register parallel backend
registerDoRedis("jobs", password = "JFEFlfki249fkjsk2~.<+JFEFl;") # the password is specific to my redis set-up
#### to stop all workers from master node:
# removeQueue("jobs")
#### from shell (to remove data on all databases)
# redis-cli -a "JFEFlfki249fkjsk2~.<+JFEFl;" flushall
#### from shell (restart redis server after flush)
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
# require(doParallel)
# registerDoParallel(cores = 1)
# registerDoParallel(cores=detectCores())
# (set.cores <- getDoParWorkers())

# --------------------------------------------------------------------------------------------
# Plot re-scaled histogram densities vs. their theoretical standard normal limit
# --------------------------------------------------------------------------------------------
make_ggplot_hist <- function(model_scen, simRes_byscen) {
  require("ggplot2")
  nsamp <- simRes_byscen[["Sim_params"]][["nsamp"]]
  Kmax <- simRes_byscen[["Sim_params"]][["Kmax"]]
  psi0 <- simRes_byscen[["psi0"]]
  # mean stats over sims:
  sims_MeanRes_byscen <- simRes_byscen[["sims_MeanRes_byscen"]]
  # raw results for all sims:
  sims_RawRes_byscen <- simRes_byscen[["sims_RawRes_byscen"]]
  # scen1_Raw_ests <- sims_RawRes_byscen[["Qh.corr"]][["est"]]
  # scen2_Raw_ests <- sims_RawRes_byscen[["Q.miss"]][["est"]]
  # scen3_Raw_ests <- sims_RawRes_byscen[["h.miss"]][["est"]]
  sigma2_0 <- sims_MeanRes_byscen[[model_scen]][,"empVar"]
  # sigma2_N <- sims_MeanRes_byscen[[model_scen]][,"ICVar"]
  all_psi_n <- sims_RawRes_byscen[[model_scen]][["est"]]

  est_names <- colnames(all_psi_n)
  all_psi_n_scaled <- (all_psi_n - psi0)
  # all_psi_n_scaled <- sqrt(nsamp/Kmax) * (all_psi_n - psi0)
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
create_OnetDAG <- function(Kmax = 5, nC = 1) {
  rcategor.int.base0 <- function(n, probs) {
    if (is.vector(probs)) {
      probs <- matrix(data = probs, nrow = n, ncol = length(probs), byrow = TRUE)
    }
    probs <- cbind(probs, 1 - rowSums(probs)) # sum each row and check if some need to be normalized
    pover1_idx <- which(probs[,ncol(probs)] < 0) # which rows of probs need to be normalized
    if (length(pover1_idx)>0) {
      warning("some categorical probabilities add up to more than 1, normalizing to add to 1")
      probs[pover1_idx, ncol(probs)] <- 0
      probs[pover1_idx, ] <- probs[pover1_idx, ] / rowSums(probs[pover1_idx, ]) # normalize
    }
    probs_cum <- matrix(nrow = nrow(probs), ncol = ncol(probs))
    probs_cum[,1] <- probs[,1]
    for (i in seq(ncol(probs))[-1]) {
      probs_cum[,i] <- rowSums(probs[,c(1:i)])
    }
    as.integer(rowSums(probs_cum - runif(nrow(probs_cum)) < 0))
  }

  `%+%` <- function(a, b) paste0(a, b)
  require("igraph")
  require("simcausal")
  options(simcausal.verbose=FALSE)
  if (nC>1) print("generating DAG with nC: " %+% nC)

  Amax <- 7
  shift <- 1
  # Categorical probabilities for exposure A conditional on value of W1:
  normprob <- function(x) x / sum(x)
  prob_Ni_W1_1 <- normprob(c(0.02, plogis(-3 - c(1:(Amax-1)) / 2)))
  prob_Ni_W1_2 <- normprob(c(0.02, plogis(-1.5 - c(1:(Amax-1)) / 3)))
  prob_Ni_W1_3 <- normprob(c(0.02, pnorm(-2*abs(2 - c(1:(Amax-1))) / 5)))
  prob_Ni_W1_4 <- normprob(c(0.02, pnorm(-2*abs(3 - c(1:(Amax-1))) / 5)))
  prob_Ni_W1_5 <- normprob(c(0.02, plogis(-4 + 2 * (c(1:(Amax-1)) - 2))))
  prob_Ni_W1_6 <- normprob(c(0.02, plogis(-4 + 2 * (c(1:(Amax-1)) - 3))))

  # --------------------------------------------------------------------------------------------
  # Network sampler(s) from igraph (regular graph model)
  # Generate regular random graphs with same degree for each node
  # Kmax - degree of each node
  # --------------------------------------------------------------------------------------------
  gen_network_k_reg <- function(n, Kmax, ...) {
    if (n < 20) Kmax <- 5
    igraph.reg <- igraph::sample_k_regular(no.of.nodes = n, k = Kmax, directed = TRUE, multiple = FALSE)
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(igraph.reg)
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat)
    return(NetInd_out$NetInd_k)
  }
  # --------------------------------------------------------------------------------------------
  # sample a network of nC independent communities, same network within each community, total sample size is n
  # --------------------------------------------------------------------------------------------
  gen_network_k_reg_indClusters <- function(n, nC = 1, Kmax, ...) {
    if (n < 20) {
      Kmax <- 5; nC <- 1
    }
    igraph.reg <- igraph::sample_k_regular(no.of.nodes = as.integer(n/nC), k = Kmax, directed = TRUE, multiple = FALSE)
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(igraph.reg)
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat)
    NetInd_mat_fixed <- NetInd_out$NetInd_k
    NetInd_mat_list <- lapply(seq(nC), function(k_comm) NetInd_mat_fixed + (k_comm - 1)*as.integer(n/nC))
    NetInd_mat <- do.call('rbind', NetInd_mat_list)
    return(NetInd_mat)
  }

  D <- DAG.empty()
  D <- D + network("NetInd_k", Kmax = Kmax, netfun = "gen_network_k_reg")

  D <- D +
      node("LatentW", distr = "rcategor.int.base0", probs = c(0.0494, 0.1823, 0.2806, 0.2680,0.1651, 0.0546))

  D <- D +
      node("W1", distr = "rcategor.int", probs = c(0.0494, 0.1823, 0.2806, 0.2680,0.1651, 0.0546)) +
      node("W2", distr = "rbern", prob = plogis(-0.2)) +
      node("WNoise", distr = "rbern", prob = plogis(-0.4)) +
      node("corrW.F1", distr = "rbern", prob = plogis(-8 + 2*LatentW + 2*LatentW[[1]])) +
      node("corrW.F2", distr = "rbern", prob = plogis(-6 + 1.5*LatentW + 1.5*LatentW[[2]])) +
      node("corrW.F3", distr = "rbern", prob = plogis(-6 + 1.5*LatentW + 1.5*LatentW[[3]])) +
      node("corrW.F4", distr = "rbern", prob = plogis(-4 + LatentW + LatentW[[4]])) +
      node("corrW.F5", distr = "rbern", prob = plogis(-4 + LatentW + LatentW[[5]]))

  # Define categorical exposure A: 1:Amax, where each A is influenced by categorical W1
  D <- D +
      node("A", distr = "rcategor.int",
          probs = (W1 == 1)*.(prob_Ni_W1_1) + (W1 == 2)*.(prob_Ni_W1_2) +
                  (W1 == 3)*.(prob_Ni_W1_3) + (W1 == 4)*.(prob_Ni_W1_4) +
                  (W1 == 5)*.(prob_Ni_W1_5) + (W1 == 6)*.(prob_Ni_W1_6)) +
      node("new.A", distr = "rconst", const = A)

  D <- D +
      node("probY", distr = "rconst",
          const = plogis(-5
                        -2*new.A +
                        # -0.10*new.A +
                        # -0.20*ifelse(nF>0, sum(new.A[[1:Kmax]])/nF, 0) +
                        # -0.5*ifelse(nF>0, sum(W2[[1:Kmax]])/nF, 0) +
                        +0.5*(W1-1) - 0.58*W2 +
                        +4*corrW.F1
                        +2*corrW.F2
                        # + 2*corrW.F3 + 2*corrW.F4 + 2*corrW.F5
                        ),
          replaceNAw0 = TRUE) +
      node("Y", distr = "rbern", prob = probY)
  D <- set.DAG(D)
}

compute_corrs <- function() {
  D <- create_OnetDAG()
  datO <- sim(D, n = 10000)
  head(datO)
  summary(datO)
  print(mean(datO$Y))
  Net_mat <- attributes(datO)$netind_cl$NetInd
  head(Net_mat)

  print("correlations b/ween W2: ")
  print(cor(datO[,"W2"], datO[,"W2"][Net_mat[,1]]))
  print("correlations b/ween corrW.Fj: ")
  print(cor(datO[,"corrW.F1"], datO[,"corrW.F1"][Net_mat[,1]]))
  print(cor(datO[,"corrW.F2"], datO[,"corrW.F2"][Net_mat[,2]]))
  print(cor(datO[,"corrW.F3"], datO[,"corrW.F3"][Net_mat[,3]]))
  print(cor(datO[,"corrW.F4"], datO[,"corrW.F4"][Net_mat[,4]]))
  print(cor(datO[,"corrW.F5"], datO[,"corrW.F5"][Net_mat[,5]]))

  print("correlations b/ween Y: ")
  print(cor(datO[,"Y"], datO[,"Y"][Net_mat[,1]]))
  print("correlations b/ween Y conditional on corrW.F1: ")
  varcorr <- "corrW.F1"
  print(cor(datO[,"Y"][datO[,varcorr]==0 & datO[,varcorr][Net_mat[,1]]==0], datO[,"Y"][Net_mat[,1]][datO[,varcorr]==0 & datO[,varcorr][Net_mat[,1]]==0]))
  print(cor(datO[,"Y"][datO[,varcorr]==1 & datO[,varcorr][Net_mat[,1]]==1], datO[,"Y"][Net_mat[,1]][datO[,varcorr]==1 & datO[,varcorr][Net_mat[,1]]==1]))

  varcorr <- "corrW.F2"
  print(cor(datO[,"Y"][datO[,varcorr]==0 & datO[,varcorr][Net_mat[,1]]==0], datO[,"Y"][Net_mat[,1]][datO[,varcorr]==0 & datO[,varcorr][Net_mat[,1]]==0]))
  print(cor(datO[,"Y"][datO[,varcorr]==1 & datO[,varcorr][Net_mat[,1]]==1], datO[,"Y"][Net_mat[,1]][datO[,varcorr]==1 & datO[,varcorr][Net_mat[,1]]==1]))

  varcorr <- "corrW.F3"
  print(cor(datO[,"Y"][datO[,varcorr]==0 & datO[,varcorr][Net_mat[,1]]==0], datO[,"Y"][Net_mat[,1]][datO[,varcorr]==0 & datO[,varcorr][Net_mat[,1]]==0]))
  print(cor(datO[,"Y"][datO[,varcorr]==1 & datO[,varcorr][Net_mat[,1]]==1], datO[,"Y"][Net_mat[,1]][datO[,varcorr]==1 & datO[,varcorr][Net_mat[,1]]==1]))

  print(cor(datO[,"Y"], datO[,"Y"][Net_mat[,2]]))
  print(cor(datO[,"Y"], datO[,"Y"][Net_mat[,3]]))
  print(cor(datO[,"Y"], datO[,"Y"][Net_mat[,4]]))
  print(cor(datO[,"Y"], datO[,"Y"][Net_mat[,5]]))


  Dset <- create_OnetDAG(Kmax = 5, nC = 1)
  shift <- 1
    Dset <- Dset +
      action("gstar",
        nodes = node("new.A", distr = "rconst", const = ifelse(A + shift <= 7, A + shift, A)),
        shift = shift)
    datFull <- sim(Dset, actions="gstar", n = 50000, rndseed = 54321)[["gstar"]]
    psi0 <- mean(datFull$Y)
    print(psi0)

}



# --------------------------------------------------------------------------------------------
# Evaluate the true value of the parameter
# (defines the action on the DAG object in Det and simulates full data)
# --------------------------------------------------------------------------------------------
eval_psi0 <- function(Kmax, shift = 1, nfull = 50000, repl = 10, nC = 1) {
  `%+%` <- function(a, b) paste0(a, b)
  Amax <- 7
  assertthat::assert_that(is.count(repl))
  if (nC>1) print("evaluting psi0 for nC: " %+% nC)
  if (repl == 1) {
    Dset <- create_OnetDAG(Kmax = Kmax, nC = nC)
    Dset <- Dset +
      action("gstar",
        nodes = node("new.A", distr = "rconst", const = ifelse(A + shift <= Amax, A + shift, A)),
        shift = shift)
    datFull <- sim(Dset, actions="gstar", n = nfull, rndseed = 54321)[["gstar"]]
    psi0 <- mean(datFull$Y)
  # parallel evaluation of psi0:
  } else { 
    ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
    ex <- c(ex, "shift", "nfull", "create_OnetDAG"); # print("ex: "); print(ex)
    psi0_repl <- foreach(t.counter = 1 : repl, .packages = c("simcausal", "igraph"), .combine = "c", .export = ex) %dopar% {
      Dset <- create_OnetDAG(Kmax = Kmax, nC = nC)
      Dset <- Dset +
        action("gstar",
          nodes = node("new.A", distr = "rconst", const = ifelse(A + shift <= Amax, A + shift, A)),
          shift = shift)
      datFull <- sim(Dset, actions="gstar", n = nfull)[["gstar"]]
      return(mean(datFull$Y))
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
create_f.gstar <- function(shift = 1) {
  Amax <- 7
  f.gstar <- function(data, ...) {
    A <- data[,"A"]
    return(ifelse(A + shift <= Amax, A + shift, A))
  }
  return(f.gstar)
}

# --------------------------------------------------------------------------------------------
# Run nsim parallel simulations for a specific scenario in Sim_params (K, nsamp, ...)
# Each run constists of 3 model specification scenarios: Q & h correct, Q misspecified, h misspecified
# Results are saved for: tmle, h_iptw (efficient IPTW), gcomp
# --------------------------------------------------------------------------------------------
runsim <- function(nsim, Sim_params, psi0, n_forpsi0, n_replpsi0) {
  # Sim_params <- All_Sims_Params[[1]]
  # --------------------------------------------------------------------------------------------
  # one run of tmlenet, output some stats
  # --------------------------------------------------------------------------------------------
  run_tmlenet <- function(datO, NetInd_mat, def_sW, def_sA, Kmax, Qform, hforms, f.gstar, psi0, bin.method, maxNperBin) {
    if (bin.method%in%"equal.mass") {
      tmlenet_options(bin.method = bin.method, maxNperBin = maxNperBin)
    } else {
      tmlenet_options(bin.method = bin.method, nbins = as.integer(nrow(datO)/maxNperBin))
    }
    # datO_input <- datO[,c("W1", "W2", "W3", "A", "Y")]
    datO_input <- datO
    hform.g0 <- hforms["hform.g0"]
    hform.gstar <- hforms["hform.gstar"]
    res <- tmlenet(data = datO_input, sW = def_sW, sA = def_sA, Anode = "A", # Ynode = "Y",
                    Kmax = Kmax, NETIDmat = NetInd_mat,
                    f_gstar1 = f.gstar, Qform = Qform,
                    hform.g0 = hform.g0,
                    hform.gstar = hform.gstar,
                    optPars = list(
                      # f_g0 = f.A, # can pass true g.0, which will then be used for sampling prior to estimating h.0
                      n_MCsims = 1)
                    )
    # psi0 estimates:
    replace_iptw_name <- rownames(res$EY_gstar1$estimates) %in% "h_iptw"

    est <- as.vector(res$EY_gstar1$estimates)
    names(est) <- rownames(res$EY_gstar1$estimates); names(est)[replace_iptw_name] <- "h.iptw" # to not crash latex
    MAD <- abs(est-psi0)
    # IC-based Var estimates:
    ICVar <- as.vector(res$EY_gstar1$vars); names(ICVar) <- rownames(res$EY_gstar1$vars); names(ICVar)[replace_iptw_name] <- "h.iptw"
    # IC-based CIs:
    CIs <- res$EY_gstar1$CIs; rownames(CIs)[replace_iptw_name] <- "h.iptw"
    cover <- (psi0 <= CIs[,2]) & (psi0 >= CIs[,1])
    CIlen <- as.vector(CIs[,2]-CIs[,1]); names(CIlen) <- rownames(CIs)
    CIlow <- as.vector(CIs[,1]); CIhigh <- as.vector(CIs[,2]);
    names(CIlow) <- names(CIhigh) <- rownames(CIs)
    return(list(est = est, MAD = MAD, cover = cover, CIlen = CIlen, CIlow = CIlow, CIhigh = CIhigh, ICVar = ICVar))
  }
  # --------------------------------------------------------------------------------------------
  # Simulation parameters
  # --------------------------------------------------------------------------------------------
  scen.name <- Sim_params[["name"]]
  nsamp <- Sim_params[["nsamp"]]
  Kmax <- Sim_params[["Kmax"]]
  shift <- Sim_params[["shift"]]
  f.gstar <- create_f.gstar(shift = shift)
  # types of models to run:
  Qforms <- Sim_params[["Qforms"]]
  hforms <- Sim_params[["hforms"]]
  bin.methods <- Sim_params[["bin.method"]]
  nobsBin_fracs <- Sim_params[["nobsBin_frac"]]
  nComm <- Sim_params[["nComm"]]
  print("running simulation with nComm: "%+%nComm)

  # loop estimation on the same simulated observed dataset changing these params:
  # **** TO DO: create this dynamically based on Qforms & hforms ****
  # Scen_params <- list(
  #   Qh.corr  = list(name = "Qh.corr", Qform = Qforms[["Q.corr"]], hforms =  hforms[["h.corr"]])
  #   # Q.miss   = list(name = "Q.miss",  Qform = Qforms[["Q.miss"]], hforms =  hforms[["h.corr"]]),
  #   # h.miss   = list(name = "h.miss",  Qform = Qforms[["Q.corr"]], hforms =  hforms[["h.miss"]])
  # )
  create_Scen_Params <- function(Q.form, Q.form.name) {
      mapply(function(bin.method, nobsBin_frac) {
                    list(name = Q.form.name%+%"_bin." %+% bin.method %+% "_nobsBinf"%+%nobsBin_frac,
                         bin.method = bin.method, nobsBin_frac = nobsBin_frac,
                         Qform = Q.form, hforms =  hforms[["h.corr"]])
                  }, rep(bin.methods, each = length(nobsBin_fracs)), rep(nobsBin_fracs, times = length(bin.methods)), SIMPLIFY = FALSE)
  }
  Scen_params_Qh.corr <- create_Scen_Params(Q.form = Qforms[["Q.corr"]], Q.form.name = "Qh.corr")
  Scen_params_Q.miss <- create_Scen_Params(Q.form = Qforms[["Q.miss"]], Q.form.name = "Q.miss")
  Scen_params <- c(Scen_params_Qh.corr, Scen_params_Q.miss)
  scen_names <- unlist(lapply(Scen_params, "[[", "name"))
  names(Scen_params) <- scen_names

  # --------------------------------------------------------------------------------------------
  # Evaluate the true value of the parameter
  # --------------------------------------------------------------------------------------------
  if (missing(psi0)) {
    message("evaluating psi0, the the true parameter value...")
    psi0 <- eval_psi0(Kmax = Kmax, shift = shift, nfull = n_forpsi0, repl = n_replpsi0, nC = nComm)
    print("psi0: " %+% psi0)
  }

  # --------------------------------------------------------------------------------------------
  # Run a parallel loop that Simulates observed data then runs tmlenet for each scenario in Scen_params
  # --------------------------------------------------------------------------------------------
  ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
  # ex <- c(ex, "rcategor.int.base0"); # print("ex: "); print(ex)
  simRes <- foreach(t.counter = 1 : nsim, .packages = c("tmlenet","simcausal", "igraph"), .export = ex) %dopar% {
    # ------------------------------------------------------------
    options(simcausal.verbose = FALSE)
    options(tmlenet.verbose = FALSE)
    # ------------------------------------------------------------
    # define the data generating distribution for the observed data:
    Dset <- create_OnetDAG(Kmax = Kmax, nC = nComm)
    datO <- sim(Dset, n = nsamp)  #, rndseed = 123456
    netind_cl <- attributes(datO)$netind_cl
    NetInd_mat <- attributes(datO)$netind_cl$NetInd

    res_byQ <- vector(mode="list", length=length(Scen_params)); names(res_byQ) <- names(Scen_params)
    for (scen.name in names(Scen_params)) {
      params <- Scen_params[[scen.name]]
      print("scen.name: " %+% scen.name); print("t.counter: "%+%t.counter)
      print(params$Qform); print(params$hforms)
      tmlenet_res <- run_tmlenet(datO = datO, NetInd_mat = NetInd_mat,
                                  def_sW = Sim_params[["def_sW"]], def_sA = Sim_params[["def_sA"]], Kmax = Kmax,
                                  Qform = params$Qform, hforms = params$hforms,
                                  f.gstar = f.gstar, psi0 = psi0,
                                  bin.method = params[["bin.method"]],
                                  maxNperBin = as.integer(params[["nobsBin_frac"]]*nsamp)
                                  )
      res_byQ[[scen.name]] <- tmlenet_res
    }
    return(res_byQ)
  }

  # --------------------------------------------------------------------------------------------
  # Collect results and take means over sims:
  # --------------------------------------------------------------------------------------------
  # to re-order results by one scenario (a list of simN results by scenario turned into a list of scenarios by simN):
  simRes_byscen <- lapply(names(Scen_params), function(scen) lapply(simRes, '[[', scen)); names(simRes_byscen) <- names(Scen_params)
  # ouput names from run_tmlenet():
  outparams <- names(simRes_byscen[[1]][[1]])
  # raw results as a list of matrices (each matrix column are results for a particular estimator):
  sims_RawRes_byscen <- lapply(simRes_byscen, function(simRes) lapply(outparams, function(outparam) t(sapply(simRes, '[[', outparam))))

  # collect mean matrix results over sims in a list (by scenario):
  sims_MeanRes_byscen <- vector(mode="list", length=length(Scen_params)); names(sims_MeanRes_byscen) <- names(Scen_params)
  for (scen in names(Scen_params)) {
    names(sims_RawRes_byscen[[scen]]) <- outparams
    # by scenario, take by column means for each output parameter returned by run_tmlenet():
    sims_MeanRes <- lapply(sims_RawRes_byscen[[scen]], base::colMeans)
    bias_abs <- abs(sims_MeanRes[["est"]] - psi0) # absolute bias
    empVar <- colVars(sims_RawRes_byscen[[scen]][["est"]]); names(empVar) <- names(sims_MeanRes[["est"]]) # empirical Var
    sims_MeanRes <- c(list(bias_abs=bias_abs, empVar=empVar), sims_MeanRes)
    sims_MeanRes <- sapply(sims_MeanRes, I)
    sims_MeanRes_byscen[[scen]] <- sims_MeanRes
  }

  sims_Summary <- list(simname = Sim_params[["name"]], psi0 = psi0, Sim_params = Sim_params,
                      sims_MeanRes_byscen = sims_MeanRes_byscen,
                      sims_RawRes_byscen = sims_RawRes_byscen
                      )
  return(sims_Summary)
}


# ********************************************************************************************
# Main function that sets-up all simulations
# ********************************************************************************************
# TO DO: Need to fix one network across all sims, evaluate psi0 also on the basis of that single network
# TO DO: Obtain the iid asymptotic var est (ignore dependence), this requires modifying tmlenet() to return the unit-specific ICs
# --------------------------------------------------------------------------------------------
run_all_sims <- function() {
  # nsim <- 5000
  nsim <- 1000
  # nsim <- 20
  # nsim <- 2
  n_forpsi0 <- 100000 # sample size for evaluting true psi.0
  n_replpsi0 <- 10
  # n_forpsi0 <- 10000 # sample size for evaluting true psi.0
  # n_replpsi0 <- 1

  # --------------------------------------------------------------------------------------------
  # Observed data sample sizes and connectivity
  # --------------------------------------------------------------------------------------------
  # simulation sample sizes:
  # nsamp_vec <- c(500, 1000, 5000, 10000, 20000, 30000, 40000, 50000)
  # nsamp_vec <- c(50000, 100000)
  # nsamp_vec <- c(25000, 50000)
  nsamp_vec <- c(500, 1000, 5000, 10000, 20000)
  # simulation network sizes (maximum number of friends for each unit):
  # Kmax_vec <- seq(5, 30, by=5)
  # Kmax_vec <- c(5)
  # Kmax_vec <- c(5, 30)
  Kmax_vec <- c(5, 10)

  # --------------------------------------------------------------------------------------------
  # Common simulation parameters (These parameters are shared across all simulations)
  # --------------------------------------------------------------------------------------------
  # Network and stochastic intervention params
  def_sW <- def.sW(W1, W2, WNoise, corrW.F1, corrW.F2, corrW.F3, corrW.F4, corrW.F5) +
            def.sW(net.W2 = ifelse(nF>0, rowSums(W2[[1:Kmax]])/nF, 0), replaceNAw0 = TRUE)
            # def.sW(W1.W2 = W1*W2) +
            # def.sW(mW1.W2 = (1-W1)*(1-W2)) +
            # def.sW(W1.W3 = W1*W3) +
            # def.sW(mW1.W3 = (1-W1)*(1-W3)) +
            # def.sW(W2.W3 = W2*W3) +
            # def.sW(mW2.W3 = (1-W2)*(1-W3)) +
            # def.sW(net.mean.W1 = ifelse(nF>0, rowSums(W1[[1:Kmax]])/nF, 0), replaceNAw0 = TRUE) +
            # def.sW(net.mean.W1.2 = ifelse(nF>0, (rowSums(W1[[1:Kmax]])/nF)^2, 0), replaceNAw0 = TRUE) +
            # def.sW(net.mean.W1.3 = ifelse(nF>0, (rowSums(W1[[1:Kmax]])/nF)^3, 0), replaceNAw0 = TRUE)

  def_sA <- def.sA(A)
            # def.sA(net.A = ifelse(nF > 0, rowSums(A[[1:Kmax]])/nF, 0), replaceNAw0 = TRUE)

  Qforms <- list(Q.corr = "Y ~ A + W1 + W2 + corrW.F1 + corrW.F2",
                 Q.miss = "Y ~ A + W1 + W2 + WNoise")

  hforms <- list(h.corr = c(hform.g0 = "A ~ W1 + W2 + corrW.F1 + corrW.F2",
                            hform.gstar = "A ~ W1 + W2 + corrW.F1 + corrW.F2"),
                 h.miss = c(hform.g0 = "A ~ W2",
                            hform.gstar = "A ~ W2"))

  # hforms <- list(h.corr = c(hform.g0 = "A + net.A ~ W1 + W2 + W3 + net.W2 + W1.eq1 + W1.eq2 + W1.eq3 + W1.eq4",
  #                           hform.gstar = "A + net.A ~ W1 + W2 + W3 + net.W2 + W1.eq1 + W1.eq2 + W1.eq3 + W1.eq4"),
  #                h.miss = c(hform.g0 = "A ~ W1 + W3",
  #                           hform.gstar = "A ~ W1 + W3"))

  Common_Sim_Params <- list(shift = 1, nComm = 1, def_sW = def_sW, def_sA = def_sA, nobsBin_frac = c(0.0025), bin.method = "equal.mass", Qforms = Qforms, hforms = hforms)
  # Common_Sim_Params <- list(shift = 1, nComm = 10, def_sW = def_sW, def_sA = def_sA, nobsBin_frac = c(0.0025, 0.005), bin.method = c("equal.mass"), Qforms = Qforms, hforms = hforms)

  # --------------------------------------------------------------------------------------------
  # Automatically populates parameter list for all combos of nsamp & Kmax in nsamp_vec & Kmax_vec
  # --------------------------------------------------------------------------------------------
  create_Sim_Params <- function(Common_Sim_Params) {
      mapply(function(nsamp, Kmax) {
                    c(
                    list(name = "N" %+% nsamp %+% ".Kmax"%+%Kmax,
                         nsamp = nsamp, Kmax = Kmax),
                    Common_Sim_Params)
                  }, rep(nsamp_vec, each = length(Kmax_vec)), rep(Kmax_vec, times = length(nsamp_vec)), SIMPLIFY = FALSE)
  }
  All_Sims_Params <- create_Sim_Params(Common_Sim_Params)
  str(All_Sims_Params)
  param_names <- lapply(All_Sims_Params, '[[', "name")
  names(All_Sims_Params) <- param_names

  # Example of a simulation parameter list:
  # All_Sims_Params <- list(
  #   N20000_Kmax5 = list(name = "N20000_Kmax5",  nsamp = 20000, Kmax = 5, shift = 1, def_sW = def_sW, def_sA = def_sA, nobsBin_frac = 0.05, bin.method = "equal.mass", Qforms = Qforms, hforms = hforms),
  #   N40000_Kmax5 = list(name = "N40000_Kmax5",  nsamp = 40000, Kmax = 5, shift = 1, def_sW = def_sW, def_sA = def_sA, nobsBin_frac = 0.05, bin.method = "equal.mass", Qforms = Qforms, hforms = hforms)
  # )

  simRes_all <- list()
  for (scenario in names(All_Sims_Params)) {
    # scenario <- "N20000_Kmax5"
    time_t <- system.time(simRes_byscen <- runsim(nsim = nsim, Sim_params = All_Sims_Params[[scenario]], n_forpsi0 = n_forpsi0, n_replpsi0 = n_replpsi0))
    print("sim time for one scenario:"); print(time_t)

    name <- simRes_byscen[["simname"]] # simulation scenario name from All_Sims_Params (Nsamp and Kmax) - should match w/ scenario var
    sims_MeanRes_byscen <- simRes_byscen[["sims_MeanRes_byscen"]] # mean result tables by model spec scenario
    sims_RawRes_byscen <- simRes_byscen[["sims_RawRes_byscen"]] # raw simulation results by model spec scenario
    model_names <- names(sims_MeanRes_byscen) # names of different model spec scenarios for Q & h:
    print("results for scenario: " %+% name); print(sims_MeanRes_byscen)
    simRes_all <- c(simRes_all, list(simRes_byscen))
    # save(list=c("simRes_byscen"), file="simRes_byscen_N40K_Kmax5.RData")

    # --------------------------------------------------------------------------------------------
    # Write sim result table to a text file
    # --------------------------------------------------------------------------------------------
    sink(file=outdir%+%"SimResults.txt", append=T)
    cat("\n#results for scenario: ", name, "\n")
    print(sims_MeanRes_byscen, quote = FALSE)
    sink()
    # --------------------------------------------------------------------------------------------
    # create ggplot histograms for each model scenario:
    # --------------------------------------------------------------------------------------------
    ggplothist_byscen <- lapply(model_names, make_ggplot_hist, simRes_byscen = simRes_byscen)
    # hist <- make_ggplot_hist(model_scen = "Qh.corr", simRes_byscen = simRes_byscen)
    names(ggplothist_byscen) <- model_names
    # plot histograms:    
    pdf(file=outdir%+%"fig_" %+% scenario %+%".pdf",width=8,height=12) # pdf
    multiplot(plotlist=ggplothist_byscen, cols=1)
    # pdf(file=outdir%+%"fig_" %+% scenario %+%".pdf",width=12,height=8) # pdf
    # multiplot(plotlist=ggplothist_byscen, cols=length(ggplothist_byscen))
    dev.off()
    # also save as tex for future latex processing:
    # require("tikzDevice")
    # tikz(file=outdir%+%"fig_" %+% scenario %+%".tex", standAlone = TRUE, width=12,height=8)
    # multiplot(plotlist=ggplothist_byscen, cols=length(ggplothist_byscen))
    # dev.off()
    # if (savetex) {
    #   system("/usr/texbin/pdflatex "%+%f_path%+%".tex")
    #   system("rm "%+%f_path%+%".tex")
    #   system("rm "%+%f_path%+%".aux")
    #   system("rm "%+%f_path%+%".log")
    # }
  }

  names(simRes_all) <- names(All_Sims_Params)
  save(list=c("simRes_all"), file=outdir%+%"simRes_all.RData")
}

run_all_sims()

