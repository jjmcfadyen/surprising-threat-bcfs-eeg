%% Preprocess EEG Data in Fieldtrip

clear all
close all
clc

dir_ft = 'D:\Toolboxes\fieldtrip-20190419\fieldtrip-20190419';
addpath(dir_ft)
ft_defaults;

dir_eeg = 'D:\Scratch\bCFS_EEG_Reanalysis\data\Exp2\EEG\Raw';
dir_behav = 'D:\Scratch\bCFS_EEG_Reanalysis\data\Exp2\Behavioural';

dir_save = 'D:\Scratch\bCFS_EEG_Reanalysis\results\eeg_preprocessing';

subjects = 1:33;

condition_labels = {'EN','UN','EF','UF'};
condition_codes = [31 32;
                   21 22
                   11 12
                   41 42];

% EEG TRIGGER CODES: 
% 1 = face onset: expected neutral
% 2 = face onset: unexpected neutral
% 3 = face onset: expected fearful
% 4 = face onset: unexpected fearful
% responses are the same as above + 20 (so: 21, 22, 23, 24)
% 99 = face offset (i.e. start of ISI)

baseline = [-.1, 0];
SL_window = [-.1 3];
RL_window = [-1 .25];

%% Prepare

% Prepare channel neighbour layout
layout = 'biosemi64.lay';

cfg = [];
cfg.method = 'template';
cfg.template = 'biosemi64_neighb.mat';
cfg.layout = layout;
cfg.feedback = 'no';

neighbours = ft_prepare_neighbours(cfg);

% Swapped channel inputs (due to broken electrodes)
swapped_channels = cell(1,33);
swapped_channels{27} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{28} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{30} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{31} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{32} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{33} = {'Nz','AFz'; 'M2','F2'};

%% Get stimulus-locked (SL) and response-locked (RL) trials

SL_trials = {};
RL_trials = {};
T = array2table(zeros(0,9), 'VariableNames', ...
    {'Subject','Trial','Condition','RT','Acc','SL_Trial','RL_Trial','SL_Artifact','RL_Artifact'});
for subj = 1:length(subjects)

    subject = subjects(subj);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else schar = ['S' num2str(subject)];
    end
    
    t = size(T,1) + 1; % start of T table for this subject

    % Load behavioural data
    load(fullfile(dir_behav,['s' num2str(subject) '_trial_info.mat'])) % log file
    
    trials = trial_info.trials';
    trials = trials(:);
    
    rt = trial_info.rt';
    rt = rt(:);
    
    acc = trial_info.acc';
    acc = acc(:);
    
    for trl = 1:length(trials)
        T.Subject{t+trl-1} = schar;
        T.Trial(t+trl-1) = trl;
        T.Condition(t+trl-1) = trials(trl,1);
        T.RT(t+trl-1) = rt(trl,1);
        T.Acc(t+trl-1) = acc(trl,1);
    end
    
    for c = 1:length(condition_labels)
        T.Condition(T.Condition == condition_codes(c,1) | T.Condition == condition_codes(c,2)) = c;
    end
    
    % Select EEG file
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
    SL_raw = ft_preprocessing(cfg);
    
    % Add trial information to T table
    if subj == 1 % started recording mid-way through trial 77...
        T.SL_Trial(1:77) = NaN;
        T.SL_Trial(78:end) = 1:length(SL_raw.trial);
    elseif length(SL_raw.trial) ~= length(trials)
        error(['Check the behavioural trials included in SL EEG file for ' schar])
    else T.SL_Trial = 1:length(SL_raw.trial);
    end
    
    %% RL Trials
    
    % Define SL trials
    cfg = [];
    cfg.dataset = filename;
    cfg.trialdef.eventtype  = 'STATUS';
    cfg.trialdef.eventvalue = [1:4]+20; % face onsets
    cfg.trialdef.prestim    = abs(RL_window(1));
    cfg.trialdef.poststim   = RL_window(end);
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

    cfg.demean = 'no'; % NO baseline correction (we'll do this manually later)
    RL_raw = ft_preprocessing(cfg);
    
    % Add trial information to T table
    if subj == 1 % started recording mid-way through trial 77...
        
        % we have to get rid of the first response because we are missing it's stimulus onset (needed for baseline correction)
        cfg = [];
        cfg.trials = 2:length(RL_raw.trial);
        RL_raw = ft_selectdata(cfg,RL_raw);
        
        T.RL_Trial(1:77) = NaN;
        
        T.RL_Trial(78:end) = -1;
        T.RL_Trial(T.RL_Trial == -1 & ~isnan(T.RT)) = 1:length(RL_raw.trial);
        T.RL_Trial(T.RL_Trial == -1) = NaN;
        
    elseif length(SL_raw.trial) ~= length(trials)
        error(['Check the behavioural trials included in RL EEG file for ' schar])
    else T.SL_Trial = 1:length(SL_raw.trial);
    end
    
    %% For both trial types...

    for e = 1:2
        
        if e == 1
            data = SL_raw;
            trial_idx = SL_raw.sampleinfo;
        elseif e == 2
            data = RL_raw;
            trial_idx = RL_raw.sampleinfo;
        end
    
        % Downsample
        cfg = [];
        cfg.resamplefs = 200;
        cfg.detrend = 'no';
        cfg.demean = 'yes';
        cfg.baselinewindow = 'all'; % use complete trial
        cfg.trials = 'all';
        cfg.sampleindex = 'no'; % whether to add a channel that has the original sample index or not
        data = ft_resampledata(cfg, data);
    
        % Correct any channel swaps
        if ~isempty(swapped_channels{subj})
            error('Need to write code to correct swapped channels')
        end

        % Interpolate noisy channels
        cfg = [];
        cfg.viewmode = 'vertical';
        ft_databrowser(cfg, data);
        
        bad_channels = input('Any bad channels to interpolate? Type {''Cz'',''Fp1''} or similar, or ENTER to skip: ');
        
        if ~isempty(bad_channels)
            cfg = [];
            cfg.badchannel = bad_channels;
            cfg.neighbours = neighbours;
            cfg.layout = layout;
            data = ft_channelrepair(cfg,raw);
        end

        % Artifact correction via ICA
        channel = ft_channelselection({'eeg','-M1','-M2'}, data);
        
        cfg = [];
        cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
        cfg.demean = 'no'; % doesn't seem to make sense for the RL epochs... (unless the demeaning is done for the entire trial?)
        cfg.channel = channel; 
        comp = ft_componentanalysis(cfg, data);
        
        cfg = [];
        cfg.component = 1:20;
        cfg.layout    = layout;
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, comp) % component topographies
        
        cfg = [];
        cfg.layout = layout;
        cfg.viewmode = 'component';
        ft_databrowser(cfg, comp) % components with waveforms
        
        reject_comp = input('Components to reject? Type [1 2 5] or ENTER to skip: ');
        
        cfg = [];
        cfg.component = reject_comp; % to be removed component(s)
        data = ft_rejectcomponent(cfg, comp, data);
        
        % Visual artifact rejection
        data.trialinfo = [1:length(data.trial)]; % change condition numbers to trial numbers so we can see which trials get rejected
        cfg = [];
        cfg.method = 'summary';
        cfg.channel = channel;
        clean = ft_rejectvisual(cfg,data);
        
        rejected = setdiff(data.trialinfo,clean.trialinfo);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % - Add the rejected list to the T table
        % - Write up the RL trials as well
        % - Make plots
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Convert to scalp current density
        cfg = [];
        cfg.method = 'spline';
        cfg.layout = layout;
        cfg.neighbours = neighbours;
        cfg.trials = 'all';
        
        scd = ft_scalpcurrentdensity(cfg,clean);
        
        if e == 1
            SL_trials{subj} = scd;
        elseif e == 2
            RL_trials{subj} = scd;
        end
        
    end
end

