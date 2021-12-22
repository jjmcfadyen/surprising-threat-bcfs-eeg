function generate_jobs_preprocess
% This script gets the template and duplicates them, according to the names
% of the 'opts' files in 'batch'
% clustertype can be 'myriad' (UCL cluster) or 'holly' (FIL cluster)

dir_home = 'D:\bCFS_EEG_Reanalysis\';
dir_batch = fullfile(dir_home,'github','scripts','cluster','batch');

addpath(dir_batch)

%% Job settings

RAM = '16G'; % RAM allocated to each job

%% Read in jobs

% Get subject labels from data folders
subjects = [1:22 24:33];
N = length(subjects);

if N==1
    tchar = '1';
else
    tchar = ['1-' num2str(N)];
end

% Read in template file
fid = fopen('template_preprocess.sh');
tline = fgetl(fid);
template = {};
while ischar(tline)
    template = [template; tline];
    tline = fgetl(fid);
end
fclose(fid);

% copy of template
F = template;

% what to call copy of template
fname = fullfile(dir_batch,['jobs_preprocess.sh']);
jobname = 'preprocess_eeg';

for i = 1:length(template) % for each line...

    L = template{i}; % get this line

    % see if there are any tags to be replaced with variables
    L = strrep(L,'[FILENAME]','preprocess_eeg');
    L = strrep(L,'[RAM]',RAM);
    L = strrep(L,'[JOBIDX]',tchar);
    L = strrep(L,'[STAGES]','[false false true true]');
    L = strrep(L,'[SUBJECTS]',['[' num2str(subjects) ']']);

    F{i} = L;

    % write file
    fid = fopen(fname,'w');
    fprintf(fid, '%s\n',F{:}) ;
    fclose(fid);
end
end
