%% Preprocess EEG Data in Fieldtrip

clear all
close all
clc

dir_ft = 'D:\Toolboxes\fieldtrip-20190419\fieldtrip-20190419';
addpath(dir_ft)
ft_defaults;

dir_eeg = 'D:\Scratch\bCFS_EEG_Reanalysis\data\Exp2\EEG\Raw';
dir_behav = 'D:\Scratch\bCFS_EEG_Reanalysis\data\Exp2\Behavioural';

subjects = 1:33;

% EEG TRIGGER CODES: 
% 1 = face onset: expected neutral
% 2 = face onset: unexpected neutral
% 3 = face onset: expected fearful
% 4 = face onset: unexpected fearful
% responses are the same as above + 20 (so: 21, 22, 23, 24)
% 99 = face offset (i.e. start of ISI)

baseline = [-.25, 0];
SL_window = [-.25 3];
RL_window = [-1 0];

%% Prepare

% Run script that has all the channel notes
badChannelList;

% Prepare channel neighbour layout
cfg = [];
cfg.method = 'template';
cfg.template = 'biosemi64_neighb.mat';
cfg.layout = 'biosemi64.lay';
cfg.feedback = 'no';

neighbours = ft_prepare_neighbours(cfg);

%% Get stimulus-locked (SL) and response-locked (RL) trials

SL_trials = {};
RL_trials = {};
for subj = 1:length(subjects)

    subject = subjects(subj);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else schar = ['S' num2str(subject)];
    end

    % Load files (log file & EEG file)
    load(fullfile(dir_behav,['s' num2str(subject) '_trial_info.mat'])) % log file
    filename = fullfile(dir_eeg,[schar '_bCFS.bdf']); % EEG file

    %% SL Trials

    % Define SL trials
    cfg = [];
    cfg.dataset = filename;
    cfg.trialdef.eventtype  = 'STATUS';
    cfg.trialdef.eventvalue = 1:4; % face onsets
    cfg.trialdef.prestim    = abs(SL_window(1));
    cfg.trialdef.poststim   = SL_window(end);
    cfg = ft_definetrial(cfg);

    % Preprocess (filter & reference & baseline)
    cfg.continuous = 'yes'; % treat it as continuous (for filtering)

    cfg.reref = 'yes';
    cfg.refchannel = 'EEG'; % 'all' = average reference

    cfg.dftfilter = 'yes';
    cfg.dftfreq = 50; % Australian line noise = 50 Hz

    cfg.bpfilter = 'yes';
    cfg.bpfreq = [.1 40]; % bandpass filter
    cfg.bpfiltord = 3; % nth order

    cfg.bsfilter = 'yes'; % bandstop filter (to remove 10 Hz mask noise)
    cfg.bsfreq = [9 11];

    cfg.demean = 'yes'; % baseline correction
    cfg.baselinewindow = baseline;
    raw = ft_preprocessing(cfg);

    % Correct any channel swaps (need to run 'badChannelList.m')
    if ~isempty(swapped_channels{subj})
        error('Need to write code to correct swapped channels')
    end

    % Interpolate noisy channels (need to run 'badChannelList.m')
    if ~isempty(bad_channels{subj})
        cfg = [];
        cfg.badchannel = bad_channels{subj};
        cfg.neighbours = neighbours;
        cfg.layout = 'biosemi64.lay';
        interp = ft_channelrepair(cfg,raw);
        data = interp;
    else data = raw;
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

    % Automatically detect EOG artifacts
    cfg = [];
    cfg.continuous = 'yes';
    cfg.artfctdef.zvalue.channel = 'VEOG';
    cfg.artfctdef.zvalue.cutoff = 4; % z-value threshold
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0.1; % in seconds
    cfg.artfctdef.zvalue.fltpadding = 0;

    cfg.artfctdef.zvalue.bpfilter = 'yes';
    cfg.artfctdef.zvalue.bpfilttype = 'but';
    cfg.artfctdef.zvalue.bpfreq = [2, 15];
    cfg.artfctdef.zvalue.bpfiltord = 4;
    cfg.artfctdef.zvalue.hilbert = 'yes';

    cfg.artfctdef.zvalue.interactive = 'yes';

    artifact_eog = ft_artifact_zvalue(cfg,data);

    % Store
    cfg = [];
    SL_trials{subj} = ft_timelockanalysis(cfg, data);
