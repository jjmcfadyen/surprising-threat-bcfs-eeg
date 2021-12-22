function preprocess_eeg(id,stages,subjects)
% This script needs to be run on a Unix machine with FSL installed
% It operates via the OSL toolbox
% Steps completed (assigned to 'stages' boolean array - e.g. [false true true false]):
%   1. Conversion (including re-referencing)
%   2. Filtering
%   3. ICA artefact removal (including downsampling to 100 Hz)
%   4. Epoching (including merging, if multiple blocks were recorded)

subject = subjects(id);
if subject < 10
    subject = ['S0' num2str(subject)];
else
    subject = ['S' num2str(subject)];
end

%% Directories

dir_osl = '/data/holly-host/jmcfadyen/osl-core-master-UNIX/';
addpath(genpath(dir_osl));
osl_startup;

dir_home = fileparts(pwd);
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
if ~exist(fullfile(dir_results,'4_epoched'))
    mkdir(fullfile(dir_results,'4_epoched'));
end

dir_output = fullfile(dir_home,'data','Exp2','eeg','tmp',subject);
if ~exist(dir_output)
    mkdir(dir_output);
end

% Channel types
eog = {'LH','RH','LV','UV'};
otherchan = {'Nz','SNz'};

% Epoch parameters
epoch_types = {'stimlocked'}; 

epoch_windows = [-0.25 3.25]; % in seconds

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

%% Convert

if stages(1)

    filelist = extractfield(dir(fullfile(dir_data,[subject '*.bdf'])),'name');
    disp(filelist)

    for r = 1:length(filelist)
    
        cd(fullfile(dir_results,'1_converted'));
        
        S             = struct;
        S.dataset     = fullfile(dir_data,filelist{r});
        S.mode        = 'continuous';
    
        D = spm_eeg_convert(S);
        
        % manual channel swaps
        swaps = {};
        if any(contains({'S27','S28','S30','S31','S32'},subject))
            swaps = {'Nz','AFz'; 'SNz','F2'};
        elseif strcmp('S33',subject)
            swaps = {'Nz','AFz'; 'M2','F2'};
        end
        if ~isempty(swaps)
            newD = D;
            for i = 1:size(swaps,1)
                
                newD = chanlabels(newD,find(strcmp(D.chanlabels,swaps{i,1})),swaps{i,2});
                newD = chanlabels(newD,find(strcmp(D.chanlabels,swaps{i,2})),swaps{i,1});
                
                newD = chantype(newD,find(strcmp(D.chanlabels,swaps{i,1})),D.chantype(find(strcmp(D.chanlabels,swaps{i,2}))));
                newD = chantype(newD,find(strcmp(D.chanlabels,swaps{i,2})),D.chantype(find(strcmp(D.chanlabels,swaps{i,1}))));
                
                newD = units(newD,find(strcmp(D.chanlabels,swaps{i,1})),D.units(find(strcmp(D.chanlabels,swaps{i,2}))));
                newD = units(newD,find(strcmp(D.chanlabels,swaps{i,2})),D.units(find(strcmp(D.chanlabels,swaps{i,1}))));
                
            end
            D = newD;
        end
        
        % edit channel types
        D = chantype(D,find(contains(D.chanlabels,eog)),'EOG');
        D = chantype(D,find(contains(D.chanlabels,otherchan)),'Other'); 

        %% Average reference
        
        if stages(2)
            S           = struct;
            S.D         = D;
            S.refchan   = 'average';
            
            D = spm_eeg_reref_eeg(S);
        end
    end
end
    
%% Filter
    
if stages(2)

    filelist = extractfield(dir(fullfile(dir_results,'1_converted',['Mspmeeg_' subject '*.mat'])),'name');
    disp(filelist)

    for r = 1:length(filelist)

        D = spm_eeg_load(fullfile(dir_results,'1_converted',filelist{r}));
        
        % notch (10 Hz)
        S           = struct;
        S.D         = D;
        S.band      = 'stop';
        S.freq      = [9.5 10.5];
        
        D = spm_eeg_filter(S);
        
        % notch (20 Hz)
        S           = struct;
        S.D         = D;
        S.band      = 'stop';
        S.freq      = [19.5 20.5];
        
        D = spm_eeg_filter(S);
        
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
    end
end

%% ICA artefact removal
    
if stages(3)

    filelist = extractfield(dir(fullfile(dir_results,'2_filtered',['fMspmeeg_' subject '*.mat'])),'name');
    disp(filelist)

    for r = 1:length(filelist)

        D = spm_eeg_load(fullfile(dir_results,'2_filtered',filelist{r}));

        % downsample to 100 Hz
        S               = [];
        S.D             = D;
        S.fsample_new   = 100; % in Hz
        D               = spm_eeg_downsample(S);

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
        delete(fullfile(dir_results,'2_filtered',['d' filelist{r}]))
    
        disp(['=================================================================='])
        disp(['=== COPY FILE SAVE LOCATION: ' fullfile(dir_results,'3_ICA',filename) ' ==='])
        disp(['=================================================================='])
    end
end

%% Epoch

if stages(4)

    filelist = extractfield(dir(fullfile(dir_results,'3_ICA',[subject '*ICA.mat'])),'name');
    disp(filelist)

    for r = 1:length(filelist)

        D = spm_eeg_load(fullfile(dir_results,'3_ICA',filelist{r}));

        % file setup
        opts           = [];
        opts.spm_files = {fullfile(dir_results,'3_ICA',filelist{r})};
        opts.datatype  = 'eeg';
        opts.dirname   = dir_output;
    
        % epoch
        opts.epoch.do    = 1;
        opts.outliers.do = 1;
    
        % switch off
        opts.coreg.do           = 0;
        opts.bad_segments.do    = 0;
        opts.downsample.do      = 0;
        opts.africa.do          = 0;
        opts.africa.todo.ica    = 0;
        opts.africa.todo.ident  = 0;
        opts.africa.todo.remove = 0;
    
        % loop for each event type in this run
        opts.epoch.trialdef = struct();
        for e = 1:length(epoch_types)
           
            switch epoch_types{e}
                case 'stimlocked'
                    for c = 1:4
                        opts.epoch.trialdef(c).conditionlabel  = [replace(epoch_types{e},'locked','') '-' condition_labels{c}]; 
                        opts.epoch.trialdef(c).eventtype        = 'STATUS';
                        opts.epoch.trialdef(c).eventvalue       = stimulus_codes(c); 
                    end
                case 'resplocked'
                    for c = 1:4
                        opts.epoch.trialdef(c).conditionlabel  = [replace(epoch_types{e},'locked','') '-' condition_labels{c}];
                        opts.epoch.trialdef(c).eventtype        = 'STATUS';
                        opts.epoch.trialdef(c).eventvalue       = response_codes(c);
                    end
            end
            
            opts.epoch.time_range                   = epoch_windows(e,:);
            
            opts = osl_run_opt(opts);
    
            disp(['=================================================================='])
            disp(['=== FILE SAVE LOCATION: ' opts.results.spm_files_epoched{1} ' ==='])
            disp(['=================================================================='])
    
            % copy file to results folder
            filename  = [subject '_r' num2str(r) '_epoched_' epoch_types{e} '.mat'];
            D = spm_eeg_load(opts.results.spm_files_epoched{1});
            D = D.copy(fullfile(dir_results,'4_epoched',filename));
    
            disp(['=================================================================='])
            disp(['=== COPY FILE SAVE LOCATION: ' fullfile(dir_results,'4_epoched',filename) ' ==='])
            disp(['=================================================================='])
            
        end
    end

    %% Merge blocks (if applicable)
%     if length(filelist) > 1
%         
%         for e = 1:length(epoch_types)
% 
%             mergelist = {};
%             for r = 1:length(filelist)
%                 mergelist{r} = fullfile(dir_results,'4_epoched',[subject '_r' num2str(r) '_epoched_' epoch_types{e} '.mat']);
%             end
%             
%             S = [];
%             S.D = mergelist;
%             merged = spm_eeg_merge(S);
%             
%             for r = 1:length(mergelist)
%                 delete(mergelist{r});
%             end
%             
%             mergename = fullfile(dir_results,'4_epoched',D.fname);
%             D = spm_eeg_load(mergename);
%             D.copy(fullfile(dir_results,'4_epoched',[subject '_r1_epoched_' epoch_types{e} '.mat']));
%             delete(mergename);
%             
%         end
%     end
end

end