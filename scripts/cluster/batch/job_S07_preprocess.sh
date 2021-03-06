#!/bin/bash -l
# Batch script to run a multi-threaded MATLAB job under SGE.
# Request X minutes of wallclock time (format hours:minutes:seconds).  --> FILLED IN BY MATLAB SCRIPT
#$ -l h_rt=2:30:00
# Request 10 gigabyte of TMPDIR space
#$ -l tmpfs=10G
# Request X gigabytes of RAM per core. --> FILLED IN BY MATLAB SCRIPT
#$ -l mem=4G
# Request a number of threads (which will use that number of cores).
# On Myriad you can set the number of threads to a maximum of 36.
#$ -pe smp 4
# Request one MATLAB licence - makes sure your job doesn't start
# running until sufficient licenses are free.
#$ -l matlab=1
# Set the name of the job.  --> FILLED IN BY MATLAB SCRIPT
#$ -N S07
# Set the working directory to somewhere in your scratch space.
# This is a necessary step as compute nodes cannot write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID.
# This directory must already exist.   --> FILLED IN BY MATLAB SCRIPT
#$ -wd /home/skgtjm6/Scratch/
# Set up directories & modules
module load xorg-utils/X11R7.7
module load matlab/full/r2018b/9.5
# These echoes output what you are about to run  --> FILLED IN BY MATLAB SCRIPT
/usr/bin/time --verbose matlab -nosplash -nodesktop -nodisplay -r "cd('/home/skgtjm6/Scratch/2021_bCFSEEG/scripts'); preprocess_eeg('S07')"
