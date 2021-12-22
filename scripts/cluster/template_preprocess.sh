#!/bin/sh

# Name for this job
#$ -N [FILENAME]

#$ -S /bin/sh

# The output from each job - the stuff you would normally see on your screen -
# is written to files. Normally one for error messages, and one for normal output.
# The line below joins the two so everything goes to one file.
#$ -j y

#$ -l vf=[RAM]
#$ -l h_vmem=[RAM]

#$ -t [JOBIDX]

time /share/apps/matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "cd('/data/holly-host/jmcfadyen/bCFS_EEG/scripts'); preprocess_eeg($SGE_TASK_ID,[STAGES],[SUBJECTS])"

