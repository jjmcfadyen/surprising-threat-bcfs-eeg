function preprocess_eeg(subject)

%% Directories

% addpath('D:\Toolboxes\spm12');
% spm('defaults','eeg')

dir_osl = '~/Scratch/2021_HrvojeProject/toolboxes/osl-core-master-UNIX/';
addpath(genpath(dir_osl));
osl_startup;

dir_home = fileparts(pwd);
addpath(fullfile(dir_home,'scripts','utils')); % contains little function scripts

dir_data = fullfile(dir_home,'data','Exp2','eeg','raw');

% Set final 'preprocessed' directory
dir_results = fullfile(dir_home,'data','Exp2','eeg','preprocessed','unfiltered');
if ~exist(dir_results)
    mkdir(dir_results)
end
if ~exist(fullfile(dir_results,'1_converted'))
    mkdir(fullfile(dir_results,'1_converted'));
end
if ~exist(fullfile(dir_results,'2_filtered'))
    mkdir(fullfile(dir_results,'2_filtered'));
end
if ~exist(fullfile(dir_results,'3_ICA'))
    mkdir(fullfile(dir_results,'3_ICA'));
end
if ~exist(fullfile(dir_results,'4_epoched_long'))
    mkdir(fullfile(dir_results,'4_epoched_long'));
end

dir_output = fullfile(dir_home,'data','Exp2','eeg','tmp',subject);
if ~exist(dir_output)
    mkdir(dir_output);
end

rawfiles = extractfield(dir(fullfile(dir_data,[subject '*.bdf'])),'name');
disp(rawfiles)

% Channel types
eog = {'LH','RH','LV','UV'};
otherchan = {'Nz','SNz'};

% Epoch parameters
epoch_types = {'fulltrial'}; % 'faceresponse' is time-locked to response (when the face can be seen), -.5 to .3
                                            % 'fulltrial' is stimulus-locked, from -.1 to 3 seconds
epoch_windows = [-3 6]; % in seconds 

% EEG TRIGGER CODES: 
% 1 = face onset: expected neutral
% 2 = face onset: unexpected neutral
% 3 = face onset: expected fearful
% 4 = face onset: unexpected fearful
% responses are the same as above + 20 (so: 21, 22, 23, 24)
% 99 = face offset (i.e. start of ISI)

condition_labels = {'EN','UN','EF','UF'};
stimulus_codes = [1 2 3 4];
response_codes = [21 22 23 24];

%% Preprocess

for r = 1:length(rawfiles)
    
%     %% Convert
%     
%     cd(fullfile(dir_results,'1_converted'));
%     
%     S             = struct;
%     S.dataset     = fullfile(dir_data,rawfiles{r});
%     S.mode        = 'continuous';
% 
%     D = spm_eeg_convert(S);
%     
%     % manual channel swaps
%     swaps = {};
%     if any(contains({'S27','S28','S30','S31','S32'},subject))
%         swaps = {'Nz','AFz'; 'SNz','F2'};
%     elseif strcmp('S33',subject)
%         swaps = {'Nz','AFz'; 'M2','F2'};
%     end
%     if ~isempty(swaps)
%         newD = D;
%         for i = 1:size(swaps,1)
%             
%             newD = chanlabels(newD,find(strcmp(D.chanlabels,swaps{i,1})),swaps{i,2});
%             newD = chanlabels(newD,find(strcmp(D.chanlabels,swaps{i,2})),swaps{i,1});
%             
%             newD = chantype(newD,find(strcmp(D.chanlabels,swaps{i,1})),D.chantype(find(strcmp(D.chanlabels,swaps{i,2}))));
%             newD = chantype(newD,find(strcmp(D.chanlabels,swaps{i,2})),D.chantype(find(strcmp(D.chanlabels,swaps{i,1}))));
%             
%             newD = units(newD,find(strcmp(D.chanlabels,swaps{i,1})),D.units(find(strcmp(D.chanlabels,swaps{i,2}))));
%             newD = units(newD,find(strcmp(D.chanlabels,swaps{i,2})),D.units(find(strcmp(D.chanlabels,swaps{i,1}))));
%             
%         end
%         D = newD;
%     end
%     
%     % edit channel types
%     D = chantype(D,find(contains(D.chanlabels,eog)),'EOG');
%     D = chantype(D,find(contains(D.chanlabels,otherchan)),'Other');
% 
%     %% Average reference
%       
%     S           = struct;
%     S.D         = D;
%     S.refchan   = 'average';
%     
%     D = spm_eeg_reref_eeg(S);
%     
    %% Filter
    
    D = spm_eeg_load(fullfile(dir_results,'1_converted',['Mspmeeg_' rawfiles{r}]));
    
%     % notch (10 Hz)
%     S           = struct;
%     S.D         = D;
%     S.band      = 'stop';
%     S.freq      = [9.5 10.5];
%     
%     D = spm_eeg_filter(S);
%     
%     % notch (20 Hz)
%     S           = struct;
%     S.D         = D;
%     S.band      = 'stop';
%     S.freq      = [19.5 20.5];
%     
%     D = spm_eeg_filter(S);
    
    % notch (50 Hz)
    S           = struct;
    S.D         = D;
    S.band      = 'stop';
    S.freq      = [49.5 50.5];
    
    D = spm_eeg_filter(S);
    
    % bandpass
    S           = struct;
    S.D         = D;
    S.band      = 'bandpass';
    S.freq      = [0.1 45];
    
    D = spm_eeg_filter(S);
    
    copy(D,fullfile(dir_results,'2_filtered',D.fname));

    %% ICA in OSL
    
    cd(dir_output);
    
    opts                                        = [];
    opts.spm_files                              = {fullfile(dir_results,'2_filtered',D.fname)};
    opts.dirname                                = fullfile(dir_results,'2_filtered');
    opts.datatype                               = 'eeg';
    
    % bad segments
    opts.bad_segments.do                        = 1;
    
    % ICA
    opts.africa.do                              = 1;
    opts.africa.todo.ica                        = 1;
    opts.africa.todo.ident                      = 1;
    opts.africa.todo.remove                     = 1;
    opts.africa.ident.max_num_artefact_comps    = 20;
    opts.africa.ident.artefact_chans            = {'EOG'};
    opts.africa.precompute_topos                = false;
    opts.africa.ident.do_kurt                   = true;
    opts.africa.ident.mains_kurt_thresh         = 0.5;
    opts.africa.ident.do_cardiac                = true;
    opts.africa.ident.do_plots                  = false;
    opts.africa.ident.do_mains                  = true;

    % switch off
    opts.downsample.do                          = 0;
    opts.highpass.do                            = 0;
    opts.epoch.do                               = 0;
    opts.outliers.do                            = 0;
    opts.coreg.do                               = 0;

    % run
    opts = osl_run_opt(opts);
    
    disp(['=================================================================='])
    disp(['=== FILE SAVE LOCATION: ' opts.results.spm_files{1} ' ==='])
    disp(['=================================================================='])

    % copy file to results folder
    filename = [subject '_r' num2str(r) '_ICA.mat'];
    D = spm_eeg_load(opts.results.spm_files{1});
    D = D.copy(fullfile(dir_results,'3_ICA',filename));

    disp(['=================================================================='])
    disp(['=== COPY FILE SAVE LOCATION: ' fullfile(dir_results,'3_ICA',filename) ' ==='])
    disp(['=================================================================='])
    
%     %% Epoch
% 
%     % file setup
%     opts           = [];
%     opts.spm_files = {fullfile(dir_results,'3_ICA',[subject '_r' num2str(r) '_ICA.mat'])};
%     opts.datatype  = 'eeg';
%     opts.dirname   = dir_output;
% 
%     % epoch
%     opts.epoch.do    = 1;
%     opts.outliers.do = 1;
% 
%     % switch off
%     opts.coreg.do           = 0;
%     opts.bad_segments.do    = 0;
%     opts.downsample.do      = 0;
%     opts.africa.do          = 0;
%     opts.africa.todo.ica    = 0;
%     opts.africa.todo.ident  = 0;
%     opts.africa.todo.remove = 0;
% 
%     % loop for each event type in this run
%     opts.epoch.trialdef = struct();
%     for e = 1:length(epoch_types)
%        
%         switch epoch_types{e}
%             case 'faceresponse'
%                 for c = 1:4
%                     opts.epoch.trialdef(c).conditionlabel  = [epoch_types{e} '-' condition_labels{c}]; % e.g. 'faceresponse-EN'
%                     opts.epoch.trialdef(c).eventtype        = 'STATUS';
%                     opts.epoch.trialdef(c).eventvalue       = response_codes(c); % matching trials (matching ev.value)
%                 end
%             case 'fulltrial'
%                 for c = 1:4
%                     opts.epoch.trialdef(c).conditionlabel  = [epoch_types{e} '-' condition_labels{c}]; % e.g. 'faceresponse-EN'
%                     opts.epoch.trialdef(c).eventtype        = 'STATUS';
%                     opts.epoch.trialdef(c).eventvalue       = stimulus_codes(c); % matching trials (matching ev.value)
%                 end
%         end
%         
%         opts.epoch.time_range                   = epoch_windows(e,:);
%         
%         opts = osl_run_opt(opts);
% 
%         disp(['=================================================================='])
%         disp(['=== FILE SAVE LOCATION: ' opts.results.spm_files_epoched{1} ' ==='])
%         disp(['=================================================================='])
% 
%         % copy file to results folder
%         filename  = [subject '_r' num2str(r) '_epochedLong_' epoch_types{e} '.mat'];
%         D = spm_eeg_load(opts.results.spm_files_epoched{1});
%         D = D.copy(fullfile(dir_results,'4_epoched_long',filename));
% 
%         disp(['=================================================================='])
%         disp(['=== COPY FILE SAVE LOCATION: ' fullfile(dir_results,'4_epoched_long',filename) ' ==='])
%         disp(['=================================================================='])
%         
%     end
end

%% Merge blocks (if applicable)
if length(rawfiles) > 1
    
    for e = 1:length(epoch_types)
        mergelist = {};
        for r = 1:length(rawfiles)
            mergelist{r} = fullfile(dir_results,'4_epoched_long',[subject '_r' num2str(r) '_epochedLong_' epoch_types{e} '.mat']);
        end
        
        S = [];
        S.D = mergelist;
        merged = spm_eeg_merge(S);
        
%         for r = 1:length(mergelist)
%             delete(mergelist{r});
%         end
%         
%         mergename = fullfile(dir_results,'4_epoched_long',D.fname);
%         D = spm_eeg_load(mergename);
%         D.copy(fullfile(dir_results,'4_epoched_long',[subject '_r1_epochedLong_' epoch_types{e} '.mat']));
%         delete(mergename);
        
    end
end

% %% Average
% 
% for e = 1:length(epoch_types)
%    
%     filename = fullfile(dir_results,'4_epoched_long',[subject '_r' num2str(r) '_epoched_' epoch_types{e} '.mat']);
% 
%     D = spm_eeg_load(filename);
%     conditionlabels = D.conditions';
%     data = ftraw(D);
%     data.conditionlabels = conditionlabels;
% 
%     save(fullfile(dir_results,'4_epoched_long_long',[subject '_avg_' epoch_types{e} '.mat']),'data');
%     
% end

%% Visualise (PC ONLY!)

% restoredefaultpath
% addpath('D:\Toolboxes\fieldtrip-20191119');
% ft_defaults;
% 
% % Look at each condition
% cmap = [0, 224, 255;
%             255, 0, 0;
%             255, 180, 29;
%             151, 34, 255]/255;
% 
% for e = 1:length(epoch_types)
% 
%     load(fullfile(dir_results,'4_epoched_long_long',[subject '_avg_' epoch_types{e} '.mat'])); % loads 'data' FT variable
%     conditionlabels = data.conditionlabels;
%     data = rmfield(data,'conditionlabels');
% 
%     thisdata = cell(1,length(condition_labels));
%     for c = 1:length(condition_labels)
%         cfg = [];
%         cfg.trials = find(contains(conditionlabels,[epoch_types{e} '-' condition_labels{c}]));
%         thisdata{c} = ft_timelockanalysis(cfg,data);
%     end
% 
%     figure
%     cfg = [];
%     cfg.layout = 'biosemi64.lay';
%     cfg.linecolor = cmap;
%     cfg.linewidth = 1.6;
%     ft_multiplotER(cfg,thisdata{:});
% 
% %         grand = ft_timelockgrandaverage([],timelock{:});
% %         figure
% %         cfg = [];
% %         cfg.layout = 'biosemi64.lay';
% %         cfg.xlim = [.1 .1];
% %         ft_topoplotER(cfg,grand);
% end

end