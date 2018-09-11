#!/bin/bash

# --------------------------------------------------------------------
# Run this script to copy simulation R files to master AWS node (host) and start the simulation
# --------------------------------------------------------------------

# change to master node ip address:
host=52.25.27.133
outdir=/home/ubuntu/sims_tmlenet
# outdir=/home/ubuntu/
# ssh ubuntu@$host

# create directory structure:
ssh ubuntu@$host 'mkdir '$outdir
ssh ubuntu@$host 'mkdir '$outdir'/Rdata_gym_boot'
# ssh ubuntu@$host 'mkdir '$outdir'/Rdata_tmp'

scp *.* ubuntu@$host:$outdir
# scp ./report/*.* ubuntu@52.11.214.104:/home/rstudio/mynewsimdir

# start a remote job:
ssh ubuntu@$host $outdir'/run_sim_script.sh&'

# after running remotely to copy back the files:
# scp ubuntu@$host:$outdir/Rdata_gym_boot/* ./Rdata_gym_boot
# scp ubuntu@52.11.214.104:$outdir/datO_input.Rdata .
# scp ubuntu@52.11.214.104:$outdir/res.Rdata .
# scp ubuntu@$host:$outdir/Rdata_catA_nC10/* ./Rdata_catA_nC10
# scp ubuntu@$host:$outdir/Rdata_catA_depW_nC1/* ./Rdata_catA_depW_nC1
# scp ubuntu@$host:$outdir/simRes_byscen.RData .



