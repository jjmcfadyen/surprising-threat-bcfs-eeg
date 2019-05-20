%% Preprocess EEG Data in Fieldtrip

clear all
close all
clc

dir_ft = 'D:\group_garrido\Jessica\fieldtrip-20190419';
addpath(dir_ft)
ft_defaults;

dir_eeg = 'G:\Private\Experiments\2_Experiment2_CFS\Final\Exp2_EEG\eeg_data_raw';
dir_behav = 'G:\Private\Experiments\2_Experiment2_CFS\Final\Exp2_EEG\behavioural_data_raw';
dir_scratch = 'D:\group_garrido\Jessica\surprising-threats-bcfs-eeg\pipeline_1';

subjects = 1:33;
conditions = {'EN','UN','EF','UF'};

% EEG TRIGGER CODES: 
% 1 = face onset: expected neutral
% 2 = face onset: unexpected neutral
% 3 = face onset: expected fearful
% 4 = face onset: unexpected fearful
% responses are the same as above + 20 (so: 21, 22, 23, 24)
% 99 = face offset (i.e. start of ISI)

baseline = [-.1, 0];
SL_window = [-.1 3];
RL_window = [-1 0.25];

RL_rt_min = abs(RL_window(1));

with_notch = true; % apply 10 Hz notch filter for mask or not (9-11 Hz)
if with_notch
    notch_tag = 'notch_';
else notch_tag = '';
end

%% Prepare

% 1st column = channel name, 2nd column = actual channel
swapped_channels = cell(1,33);
swapped_channels{27} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{28} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{30} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{31} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{32} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{33} = {'Nz','AFz'; 'M2','F2'};

% Prepare channel neighbour layout
layout = 'biosemi64.lay';

cfg = [];
cfg.method = 'template';
cfg.template = 'biosemi64_neighb.mat';
cfg.layout = layout;
cfg.feedback = 'no';

neighbours = ft_prepare_neighbours(cfg);

%% Get stimulus-locked (SL) and response-locked (RL) trials

SL_trials = {};
RL_trials = {};
SL_rejected = {};
RL_rejected = {};
for subj = 1:length(subjects)

    subject = subjects(subj);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else schar = ['S' num2str(subject)];
    end

    % Load files (log file & EEG file)
    load(fullfile(dir_behav,['s' num2str(subject) '_trial_info.mat'])) % log file
    filename = fullfile(dir_eeg,[schar '_bCFS.bdf']); % EEG file

    %% Preprocess

    % Define SL trials
    cfg = [];
    cfg.dataset = filename;
    cfg.trialdef.eventtype  = 'STATUS';
    cfg.trialdef.eventvalue = 1:4; % face onsets
    cfg.trialdef.prestim    = abs(SL_window(1));
    cfg.trialdef.poststim   = SL_window(end);
    cfg = ft_definetrial(cfg);
    
    cfg.dataset = filename;
    cfg.continuous = 'yes'; % treat it as continuous (for filtering)

    cfg.reref = 'yes';
    cfg.refchannel = 'EEG'; % 'all' = average reference

    cfg.dftfilter = 'yes';
    cfg.dftfreq = 50; % Australian line noise = 50 Hz

    cfg.bpfilter = 'yes';
    cfg.bpfreq = [.1 40]; % bandpass filter
    cfg.bpfiltord = 3; % nth order (3 is the only one that worked for S01)

    if with_notch
        cfg.bsfilter = 'yes'; % bandstop filter (to remove 10 Hz mask noise)
        cfg.bsfreq = [9 11];
    end

    cfg.demean = 'yes'; % baseline correction
    cfg.baselinewindow = baseline;
    SL_raw = ft_preprocessing(cfg);
    
    % Define RL trials
    cfg = [];
    cfg.dataset = filename;
    cfg.trialdef.eventtype  = 'STATUS';
    cfg.trialdef.eventvalue = 21:24; % response onsets
    cfg.trialdef.prestim    = abs(RL_window(1));
    cfg.trialdef.poststim   = RL_window(end);
    cfg = ft_definetrial(cfg);
    
    cfg.dataset = filename;
    cfg.continuous = 'yes'; % treat it as continuous (for filtering)

    cfg.reref = 'yes';
    cfg.refchannel = 'EEG'; % 'all' = average reference

    cfg.dftfilter = 'yes';
    cfg.dftfreq = 50; % Australian line noise = 50 Hz

    cfg.bpfilter = 'yes';
    cfg.bpfreq = [.1 40]; % bandpass filter
    cfg.bpfiltord = 3; % nth order (3 is the only one that worked for S01)

    if with_notch
        cfg.bsfilter = 'yes'; % bandstop filter (to remove 10 Hz mask noise)
        cfg.bsfreq = [9 11];
    end

    cfg.demean = 'no'; % baseline correction
    RL_raw = ft_preprocessing(cfg);
    
    % %% Use SL baseline for RL
    % first need to find the SL trials where a response was made...
    trials = trial_info.trials';
    trials = trials(:);
    rt = trial_info.rt';
    rt = rt(:);
    
    if subj == 1

        trial_idx = 78:length(trials);
        trials = trials(trial_idx);
        rt = rt(trial_idx);
        
        % delete the first response trial, as we can't put a baseline to it
        cfg = [];
        cfg.trials = 2:length(RL_raw.trial);
        RL_raw = ft_selectdata(cfg,RL_raw);
        
    else trial_idx = 1:length(trials);
    end

    response_idx = find(~isnan(rt)); % this should equal the no. of triggers in RL_raw
    
    for trl = 1:length(response_idx)
        stim_trial = response_idx(trl);
        baseline_idx = 1:find(SL_raw.time{stim_trial} == 0);
        this_baseline = mean(SL_raw.trial{stim_trial}(:,baseline_idx),2); % averaged over baseline time window
        RL_raw.trial{trl} = RL_raw.trial{trl} - this_baseline;
    end
    
    % Put trials into easy to understand table
    response_idx = nan(length(rt),1);
    response_idx(~isnan(rt),1) = 1:sum(~isnan(rt));
    
    trial_summary = array2table([trial_idx',rt,SL_raw.trialinfo,[1:length(SL_raw.trial)]',response_idx],...
        'VariableNames',{'Behavioural_Trials','RT','Condition','SL_Trials','RL_Trials'});    
    
    trial_summary.RL_Rejected = zeros(size(trial_summary,1),1);
    trial_summary.RL_Rejected(isnan(response_idx)) = nan;
    trial_summary.SL_Rejected = zeros(size(trial_summary,1),1);
    for e = 1:2 % for both the SL and RL trials...
        
        if e == 1
            data = SL_raw;
            epoch_tag = 'SL';
        elseif e == 2
            data = RL_raw;
            epoch_tag = 'RL';
        end
        
        % Correct any channel swaps
        if ~isempty(swapped_channels{subj})
            error('Need to write code to correct swapped channels')
        end

        % Extract VEOG information, re-reference, and add it back into the dataset
        cfg = [];
        cfg.channel = {'UV','LV'}; % get the upper and lower vertical EOG channels
        cfg.reref = 'yes';
        cfg.refchannel = {'UV'}; % reference them to each other (doesn't matter which one)
        veog = ft_preprocessing(cfg,data);

        cfg = [];
        cfg.channel = 'LV'; % pick one of them
        veog = ft_selectdata(cfg,veog);
        veog.label = {'VEOG'}; % relabel

        cfg = [];
        cfg.channel = 1:64; % only select the EEG channels
        data = ft_selectdata(cfg,data);

        cfg = [];
        data = ft_appenddata(cfg, data, veog);

        raw = data;
        save(fullfile(dir_scratch,[schar '_' epoch_tag '_' notch_tag 'raw.mat']),'raw');

        % Visually reject noisy channels (loop through trials)
        cfg = [];
        cfg.viewmode = 'vertical';
        ft_databrowser(cfg,data); % visually check

        bad_chan = input('Channels to interpolate? (type {''-Cpz'',''-F1''} OR press ENTER to skip): ');

%         cfg = []; % INTERPOLATE
%         cfg.channel = bad_chan;
%         data = ft_selectdata(cfg,data);
        
        % Use ICA to remove EOG artifacts (and others)
        cfg = [];
        cfg.method = 'runica';
        cfg.channels = 1:64;
        comp = ft_componentanalysis(cfg,data); % TAKES 3-4 HOURS FOR SL AND 1 HOUR FOR RL
        save(fullfile(dir_scratch,[schar '_' epoch_tag '_' notch_tag 'comp.mat']),'comp');

        cfg = [];
        cfg.component = 1:20;
        cfg.layout = layout;
        cfg.comment = 'no';
        ft_topoplotIC(cfg,comp); % view topographies

        cfg = [];
        cfg.layout = layout;
        cfg.viewmode = 'component';
        ft_databrowser(cfg, comp) % view waveforms as well

        art_comp = input('Components to delete? (type [1 2] OR press ENTER to skip): ');

        cfg = [];
        cfg.component = art_comp;
        clean = ft_rejectcomponent(cfg, comp, data); % remove from data

        % Visual artifact rejection (amplitude, variance)
        cfg = [];
        cfg.method = 'summary';
        cfg.channel = 1:64;
        cfg.neighbours = neighbours;
        cfg.layout = layout;
        cfg.keeptrial = 'no';
        clean = ft_rejectvisual(cfg,clean);

        save(fullfile(dir_scratch,[schar '_' epoch_tag '_' notch_tag 'clean.mat']),'clean');
        
        % Get index of rejected trials
        artfctdef = clean.cfg.artfctdef.summary.artifact;
        if e == 1
            idx = find(~isnan(trial_summary.SL_Trials));
            for i = 1:size(artfctdef,1)
                trial_summary.SL_Rejected(idx(find(data.sampleinfo == artfctdef(i,1))),1) = 1;
            end
        elseif e == 2
            idx = find(~isnan(trial_summary.RL_Trials));
            for i = 1:size(artfctdef,1)
                trial_summary.RL_Rejected(idx(find(data.sampleinfo == artfctdef(i,1))),1) = 1;
            end
        end
        
        % Convert to scalp current density
        cfg = [];
        cfg.method = 'spline';
        cfg.layout = layout;
        cfg.neighbours = neighbours;
        scd = ft_scalpcurrentdensity(cfg, clean);
        save(fullfile(dir_scratch,[schar '_' epoch_tag '_' notch_tag 'scd.mat']),'scd');
        
        % Get data per condition
        cfg = [];
        cfg.channel = 1:64;
        cdata = {};
        for c = 1:length(conditions)
            if e == 1
                cfg.trials = find(scd.trialinfo == c);
            elseif e == 2
                this_rt_idx = trial_summary.RT(trial_summary.RL_Rejected == 0);
                cfg.trials = find(scd.trialinfo == c+20 & this_rt_idx >= RL_rt_min);
                disp([num2str(length(cfg.trials)) ' trials identified for condition ' num2str(c) ' with RT >= ' num2str(RL_rt_min) ' seconds'])
            end
            cdata{c} = ft_timelockanalysis(cfg,scd);
        end
        
        % View condition averages
        cfg = [];
        cfg.layout = layout;
        cfg.interactive = 'yes';
        cfg.showoutline = 'yes';
        ft_multiplotER(cfg, cdata{:});
        
        if e == 1
            for c = 1:length(conditions)
                SL_trials{subj,c} = cdata{c};
            end
            save(fullfile(dir_scratch,['SL_' notch_tag 'trials.mat']),'SL_trials');
        elseif e == 2
            for c = 1:length(conditions)
                RL_trials{subj,c} = cdata{c};
            end
            save(fullfile(dir_scratch,['RL_' notch_tag 'trials.mat']),'RL_trials');
        end

    end
    
    save(fullfile(dir_scratch,[schar '_' notch_tag 'trial_summary.mat'),'trial_summary');
    
    %% 
    
    
    
end
