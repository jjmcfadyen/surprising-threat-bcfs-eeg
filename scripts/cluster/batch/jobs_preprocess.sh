#!/bin/sh
# Name for this job
#$ -N preprocess_eeg
#$ -S /bin/sh
# The output from each job - the stuff you would normally see on your screen -
# is written to files. Normally one for error messages, and one for normal output.
# The line below joins the two so everything goes to one file.
#$ -j y
#$ -l vf=16G
#$ -l h_vmem=16G
#$ -t 1-32
time /share/apps/matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "cd('/data/holly-host/jmcfadyen/bCFS_EEG/scripts'); preprocess_eeg($SGE_TASK_ID,[false false true true],[1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  24  25  26  27  28  29  30  31  32  33])"
