#!/bin/bash
export R_VSIZE=1GB
cd /home/ubuntu/sims_tmlenet
Rscript --no-save --no-restore --verbose --vanilla --slave simulation_gym_depW.R > simulation_gym_depW.Rout 2>&1
echo "script execution complete"
exit 0