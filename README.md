# Surprising Threats in bCFS: Behavioural & EEG Experiments
This repository contains the experiment and analysis scripts for our manuscript, "Surprising threats accelerate evidence accumulation for conscious perception," by Jessica McFadyen, Cooper Smout, Naotsugu Tsuchiya, Jason B. Mattingley, & Marta I. Garrido. It is currently available at bioRxiv: https://doi.org/10.1101/525519

McFadyen, J., Smout, C., Tsuchiya, N., Mattingley, J. B., & Garrido, M. I. (2019). Surprising threats accelerate evidence accumulation for conscious perception. bioRxiv, 525519.

## The Experiment
There were two experiments in this study: one that was behavioural only, and another that incorporated EEG recordings. Both experiments used a breaking continuous flash suppression (bCFS) paradigm, where face stimuli were unconsciously presented via dichoptic presentation (i.e. a flashing mask in one eye, a face in the other eye). Participants had to report the orientation of the face image (i.e. tilted 5 degrees either to the left or the right) as quickly and accurately as possible. To do this, they had to wait until the face became visible to them - hence, response times should indicate the earliest reportable time of conscious discrimination.

The experimental paradigm requires PsychToolbox version 3 and a set of face stimuli. We used face stimuli from the following databases and edited them according to the procedure outlined in the manuscript:
* Karolinska Directed Emotional Faces (KDEF): http://www.emotionlab.se/kdef/
* Amsterdam Dynamic Facial Expressions Set (ADFES): https://aice.uva.nl/research-tools/adfes-stimulus-set/download-login-required/download.html
* NimStim: original website referred to in Tottenham et al. 2009 paper no longer works, so could use this instead - https://danlab7.wixsite.com/nimstim
* Warsaw Set of Emotional Facial Expression Pictures (WSEFEP): http://www.emotional-face.org

## The Data
All raw and processed data from both experiments can be found on the accompanying Open Science Framework project page: https://osf.io/p3du5/

Behavioural data are stored as MATLAB .mat files and EEG data are stored as .BDF files. All data does not contain any identifying information.

## The Analysis

