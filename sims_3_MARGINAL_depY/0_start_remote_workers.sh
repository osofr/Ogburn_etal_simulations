#!/bin/bash

# -------------------------------------------------------------------------------------------------------------------
# Run this script to automatically initiate redis client (R worker jobs) for all AWS nodes
# -------------------------------------------------------------------------------------------------------------------

# change to the ip address of the master node:
master=52.25.27.133

robins=52.37.230.214
monikas=52.37.117.22
# when using input arg to spec the ip of the master node:
inputip=$1

# -------------------------------------------------------------------------------------------------------------------
# NOTE: Requires first installing AWS Command Line Interface at:
# http://aws.amazon.com/cli/
# -------------------------------------------------------------------------------------------------------------------
# describes all active instances:
aws ec2 describe-instances --query 'Reservations[*].Instances[*].[Placement.AvailabilityZone, State.Name, InstanceType, InstanceId, PublicIpAddress, PublicDnsName]' --output table
# aws ec2 describe-instances --query 'Reservations[*].Instances[*].[PublicDnsName]' --output text
# aws ec2 describe-instances --query 'Reservations[*].Instances[*].[PublicIpAddress]' --output text

if [[ -n "$inputip" ]]; then
  master=$inputip
  echo "input master ip:"
  echo "$master"
else
  echo "default master ip:"
  echo "$master"
fi

hosts=$(aws ec2 describe-instances --output text --query 'Reservations[*].Instances[*].[PublicIpAddress]')

for host in $hosts; do
  echo "setting up host:"$host
  if [ $host == $robins ] || [ $host == $monikas ]
  then
    echo "robins host $robins, do nothing"
  else
    if [ $host == $master ]; then
      echo "host is equal to master $master, setting up remote"
      # copy init R script to all available nodes:
      scp -o StrictHostKeyChecking=no ./_init_local_workers.R ubuntu@$host:~/
      # remotely call R script to init workers:
      ssh ubuntu@$host 'Rscript --no-save --no-restore --verbose --vanilla --slave ~/_init_local_workers.R' $master '> ~/_init_local_workers.Rout 2>&1 &'
      sleep 0.5
    fi
    if [ $host != $master ]; then
      echo "host $host is not equal to master $master"
      # copy init R script to all available nodes:
      scp -o StrictHostKeyChecking=no ./_init_local_workers.R ubuntu@$host:~/
      # remotely call R script to init workers:
      ssh ubuntu@$host 'Rscript --no-save --no-restore --verbose --vanilla --slave ~/_init_local_workers.R' $master '> ~/_init_local_workers.Rout 2>&1 &'
      # ssh ubuntu@52.26.37.64 'Rscript --no-save --no-restore --verbose --vanilla --slave ~/_init_local_workers.R' $master '> ~/_init_local_workers.Rout 2>&1 &'
      sleep 0.5
    fi
  fi
done