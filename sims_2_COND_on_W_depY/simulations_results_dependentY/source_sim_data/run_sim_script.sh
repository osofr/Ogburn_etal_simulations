#!/bin/bash
Rscript --no-save --no-restore --verbose --vanilla --slave simulation_gym_depY_v3.R > simulation_gym_depY_v3.Rout 2>&1
echo "script execution complete"
exit 0

