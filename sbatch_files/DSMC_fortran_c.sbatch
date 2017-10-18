#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name  
#SBATCH --job-name=DSMC

#################  
#a file for job output, you can check job progress, append the job ID with %j to make it unique
#SBATCH --output=out/DSMC.%j.out

#################
# a file for errors from the job
#SBATCH --error=err/DSMC.%j.err

#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#max is 2 days
#SBATCH --time=4:00:00


#################
#Quality of Service (QOS); ~= job priority
#   --qos=long -  max job length of 7 days
#   --qos=normal -  up to 48 hours 
#SBATCH --qos=normal

#################
# choose partition to submit to (dev, normal, gpu, owners, hns, bigmem (jobs requiring >64GB RAM)
#SBATCH -p normal 

#################
#number of nodes you are requesting
#SBATCH --nodes=1

################
#number of processors?

#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=4000

#################
# Have SLURM send you an email when the job ends or fails
#SBATCH --mail-type=END,FAIL # notifications for job done & fail 
#SBATCH --mail-user=wriedel@stanford.edu

#-----------------------------------------------------------------------------
#now run normal batch commands
# note the "CMD BATCH is an R specific command
module load python/2.7.5

cd /home/wriedel/DSMC_FORTRAN

cp input_parameter_files/Input_Parameters_full_1D_c Input/Input_Parameters

rm Main
make
./Main








