#!/usr/bin/Rscript --vanilla
# R CMD BATCH --no-save Rset_up.R Rset_up.Rout

# ====================================================================================
# This script is used to initialize R redis on a remote worker server
# Can be also used to install any specific R packages
# ====================================================================================

# ------------------------------------------------------------------------------------
# controlling redis server on local MAC OS X:
# ------------------------------------------------------------------------------------
# To have launchd start redis at login:
#   ln -sfv /usr/local/opt/redis/*.plist ~/Library/LaunchAgents
# Then to load redis now:
#   launchctl load ~/Library/LaunchAgents/homebrew.mxcl.redis.plist
# Or, if you don't want/need launchctl, you can just run:
#   redis-server /usr/local/etc/redis.conf
# to test the server:
# redis-cli ping
# then in R to start one worker:
# library("doRedis");
# library("parallel");
# startLocalWorkers(n=1,queue="jobs")
# ------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
print("received args: "); print(args)

if (length(args)==0) stop("need ip address for the master node")
master <- args[1]

options(repos = "http://cran.cnr.berkeley.edu/")
print(.libPaths())
devtools::install_github('osofr/simcausal', ref = "master", build_vignettes = FALSE)
devtools::install_github('osofr/tmlenet', ref = "master", build_vignettes = FALSE)
# install.packages("tikzDevice")
# install.packages("locfit")
# install.packages("bigmemory")
# install.packages("biganalytics")
# install.packages("biganalytics")
# install.packages("/home/ubuntu/tmlenet.deprecated_0.0.0.9000.tar.gz", repos = NULL, type="source", dependencies=TRUE)

# ------------------------------------------------------------------------------------
# WORKER SET-UP SCRIPT (needs to be executed on each worker node)
# The ip address is for the single master node with redis server running
# (the master node is the one on which the simulation R code is to be exectuted)
# ------------------------------------------------------------------------------------
library("doRedis");
library("parallel");
startLocalWorkers(n=detectCores(),queue="jobs",host=master, password="JFEFlfki249fkjsk2~.<+JFEFl;")

