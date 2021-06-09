function generate_jobs_preprocess
% This script gets the template and duplicates them, according to the names
% of the 'opts' files in 'batch'
% clustertype can be 'myriad' (UCL cluster) or 'holly' (FIL cluster)

dir_home = 'D:\bCFS_EEG_Reanalysis\';
dir_batch = fullfile(dir_home,'scripts','cluster','batch');

addpath(dir_batch)

%% Job settings

timechar = '2:30:00'; % max job duration (hours : minutes : seconds)
RAM = '4G'; % RAM allocated to each job

%% Read in jobs

% Get subject labels from data folders
subjects = [1:22 24:33];

N = length(subjects);

% Read in template file
fid = fopen('template_preprocess.sh');
tline = fgetl(fid);
template = {};
while ischar(tline)
    template = [template; tline];
    tline = fgetl(fid);
end
fclose(fid);

% For each subject, determine number of runs and duplicate the template
for s = 1:N
    
    if subjects(s) < 10
        subject = ['S0' num2str(subjects(s))];
    else
        subject = ['S' num2str(subjects(s))];
    end

    % copy of template
    F = template;

    % what to call copy of template
    fname = fullfile(dir_batch,['job_' subject '_preprocess.sh']);
    jobname = [subject '_preprocess'];

    for i = 1:length(template) % for each line...

        L = template{i}; % get this line

        % see if there are any tags to be replaced with variables
        L = strrep(L,'[TIME]',timechar);
        L = strrep(L,'[RAM]',RAM);
        L = strrep(L,'[SUBJECT]',subject);

        F{i} = L;

        % write file
        fid = fopen(fname,'w');
        fprintf(fid, '%s\n',F{:}) ;
        fclose(fid);
    end
end
end
