%% Preprocess EEG Data in Fieldtrip

clear all
close all
clc

dir_ft = 'D:\group_garrido\Jessica\fieldtrip-20190419';
addpath(dir_ft)
ft_defaults;

dir_eeg = 'D:\group_garrido\Jessica\bcfs_eeg_reanalysis\data\Exp2\EEG';
dir_behav = 'D:\group_garrido\Jessica\bcfs_eeg_reanalysis\data\Exp2\Behavioural';

dir_save = 'D:\group_garrido\Jessica\bcfs_eeg_reanalysis\results';

subjects = [3 26]; %setdiff(1:33,[3 23 26]);

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

if ~exist(fullfile(dir_save,'T_table.mat'))
    T = array2table(zeros(0,9), 'VariableNames', ...
    {'Subject','Trial','Condition','RT','Acc','SL_Trial','RL_Trial','SL_Artifact','RL_Artifact'});
else load(fullfile(dir_save,'T_table.mat'));
end

for subj = 1:length(subjects)

    close all
    
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
        T.Subject(t+trl-1) = subject;
        T.Trial(t+trl-1) = trl;
        T.Condition(t+trl-1) = trials(trl,1);
        T.RT(t+trl-1) = rt(trl,1);
        T.Acc(t+trl-1) = acc(trl,1);
    end
    
    for c = 1:length(condition_labels)
        T.Condition(T.Condition == condition_codes(c,1) | T.Condition == condition_codes(c,2)) = c;
    end
    
    % Select EEG file
    if subject == 26
        filename = {fullfile(dir_eeg,[schar '_bCFS_b1.bdf']),fullfile(dir_eeg,[schar '_bCFS_b2.bdf'])};
    else
        filename = fullfile(dir_eeg,[schar '_bCFS.bdf']); % EEG file
    end

    %% SL Trials

    if subject == 26
        bdata = {};
        for b = 1:2
            cfg = [];
            cfg.dataset = filename{b};
            cfg.trialdef.eventtype  = 'STATUS';
            cfg.trialdef.eventvalue = 1:4; % face onsets
            cfg.trialdef.prestim    = abs(SL_window(1));
            cfg.trialdef.poststim   = SL_window(end);
            cfg = ft_definetrial(cfg);
            bdata{b} = ft_preprocessing(cfg);
        end
        cfg = [];
        SL_raw = ft_appenddata(cfg,bdata{1},bdata{2});
    else
    
        % Define SL trials
        cfg = [];
        cfg.dataset = filename;
        cfg.trialdef.eventtype  = 'STATUS';
        cfg.trialdef.eventvalue = 1:4; % face onsets
        cfg.trialdef.prestim    = abs(SL_window(1));
        cfg.trialdef.poststim   = SL_window(end);
        cfg = ft_definetrial(cfg);
        
    end

    % Preprocess (filter & reference & baseline)
    cfg.continuous = 'yes'; % treat it as continuous (for filtering)

    cfg.reref = 'yes';
    cfg.refchannel = 'EEG'; % 'all' = average reference

    cfg.dftfilter = 'yes';
    cfg.dftfreq = 50; % Australian line noise = 50 Hz

    cfg.bpfilter = 'yes';
    cfg.bpfreq = [.5 30]; % bandpass filter
    cfg.bpfiltord = 3; % nth order

    cfg.bsfilter = 'yes'; % bandstop filter (to remove 10 Hz mask noise)
    cfg.bsfreq = [9 11];

    cfg.demean = 'yes'; % baseline correction
    cfg.baselinewindow = baseline;
    
    if subject == 26
        SL_raw = ft_preprocessing(cfg,SL_raw);
    else
        SL_raw = ft_preprocessing(cfg);
    end
    
    % add the 20 Hz harmonic notch filter
    cfg = [];
    cfg.bsfilter = 'yes'; % bandstop filter (to remove 10 Hz mask noise)
    cfg.bsfreq = [19 21];
    SL_raw = ft_preprocessing(cfg,SL_raw);
    
    % Add trial information to T table
    if subject == 1 % started recording mid-way through trial 77...
        T.SL_Trial(T.Subject == 11 & T.Trial < 78) = NaN;
        T.SL_Trial(T.Subject == 11 & T.Trial > 77) = 1:length(SL_raw.trial);
    elseif subject == 18 % trigger for last EEG trial is missing
        T.SL_Trial(T.Subject == 18 & T.Trial == 1260) = NaN;
        T.SL_Trial(T.Subject == 18 & ~isnan(T.SL_Trial)) = 1:length(SL_raw.trial);
    elseif length(SL_raw.trial) ~= length(trials)
        error(['Check the behavioural trials included in SL EEG file for ' schar])
    else
        T.SL_Trial(t:end) = 1:length(SL_raw.trial);
    end
    
    %% RL Trials
    
    if subject == 26
        bdata = {};
        for b = 1:2
            cfg = [];
            cfg.dataset = filename{b};
            cfg.trialdef.eventtype  = 'STATUS';
            cfg.trialdef.eventvalue = [[1:4]+20 [1:4]+30]; % face onsets
            cfg.trialdef.prestim    = abs(RL_window(1));
            cfg.trialdef.poststim   = RL_window(end);
            cfg = ft_definetrial(cfg);
            bdata{b} = ft_preprocessing(cfg);
        end
        cfg = [];
        RL_raw = ft_appenddata(cfg,bdata{1},bdata{2});
    else
        % Define RL trials
        cfg = [];
        cfg.dataset = filename;
        cfg.trialdef.eventtype  = 'STATUS';
        cfg.trialdef.eventvalue = [[1:4]+20 [1:4]+30]; % face onsets
        cfg.trialdef.prestim    = abs(RL_window(1));
        cfg.trialdef.poststim   = RL_window(end);
        cfg = ft_definetrial(cfg);
    end

    % Preprocess (filter & reference & baseline)
    cfg.continuous = 'yes'; % treat it as continuous (for filtering)

    cfg.reref = 'yes';
    cfg.refchannel = 'EEG'; % 'all' = average reference

    cfg.dftfilter = 'yes';
    cfg.dftfreq = 50; % Australian line noise = 50 Hz

    cfg.bpfilter = 'yes';
    cfg.bpfreq = [.5 30]; % bandpass filter
    cfg.bpfiltord = 3; % nth order

    cfg.bsfilter = 'yes'; % bandstop filter (to remove 10 Hz mask noise)
    cfg.bsfreq = [9 11];

    cfg.demean = 'no'; % NO baseline correction (we'll do this manually later)
    
    if subject == 26
        RL_raw = ft_preprocessing(cfg,RL_raw);
    else
        RL_raw = ft_preprocessing(cfg);
    end
        
    % add the 20 Hz harmonic notch filter
    cfg = [];
    cfg.bsfilter = 'yes'; % bandstop filter (to remove 10 Hz mask noise)
    cfg.bsfreq = [19 21];
    RL_raw = ft_preprocessing(cfg,RL_raw);
    
    % Add trial information to T table
    if subject == 1 % started recording mid-way through trial 77...
        
        % we have to get rid of the first response because we are missing it's stimulus onset (needed for baseline correction)
        cfg = [];
        cfg.trials = 2:length(RL_raw.trial);
        RL_raw = ft_selectdata(cfg,RL_raw);
        
        T.RL_Trial(1:77) = NaN;
        
        T.RL_Trial(78:end) = -1;
        T.RL_Trial(T.RL_Trial == -1 & ~isnan(T.RT)) = 1:length(RL_raw.trial);
        T.RL_Trial(T.RL_Trial == -1) = NaN;
        
    elseif subject == 18
        
        % we have to get rid of the last response because we are missing it's stimulus onset (needed for baseline correction)
        cfg = [];
        cfg.trials = 1:length(RL_raw.trial)-1;
        RL_raw = ft_selectdata(cfg,RL_raw);
        
        % get onsets of responses to match
        for trl = 1:size(RL_raw.sampleinfo,1)
            this_rt_samp = RL_raw.sampleinfo(trl,2) - RL_window(2)*RL_raw.fsample;
            tmp = SL_raw.sampleinfo(:,2) - this_rt_samp;
            idx = find(tmp > 0 & tmp < 3*RL_raw.fsample);
            if length(idx) > 1
                error(['Something went wrong at line 226 for ' schar])
            end
            T.RL_Trial(t+idx-1) = trl;
        end
        T.RL_Trial(T.RL_Trial == 0) = NaN;
        
        
    elseif length(SL_raw.trial) ~= length(trials)
        error(['Check the behavioural trials included in SL EEG file for ' schar])
    elseif subject == 26
        T.RL_Trial(T.Subject == subject & ~isnan(T.RT)) = 1:length(RL_raw.trial);
    else
        
        % get onsets of responses to match
        for trl = 1:size(RL_raw.sampleinfo,1)
            this_rt_samp = RL_raw.sampleinfo(trl,2) - RL_window(2)*RL_raw.fsample;
            tmp = SL_raw.sampleinfo(:,2) - this_rt_samp;
            idx = find(tmp > 0 & tmp < 3*RL_raw.fsample);
            if length(idx) > 1
                error(['Something went wrong at line 226 for ' schar])
            end
            T.RL_Trial(t+idx-1) = trl;
        end
        T.RL_Trial(T.RL_Trial == 0) = NaN;
        
    end
    
    % Baseline correct RL trials using corresponding SL pre-stimulus period
    rl_trials = find(~isnan(T.RL_Trial(t:end)));
    for trl = 1:length(rl_trials)
        corresponding_SL = T.SL_Trial(t+rl_trials(trl)-1);
        tw = 1:find(SL_raw.time{corresponding_SL} == 0);% in samples
        RL_raw.trial{trl}(1:64,:) = RL_raw.trial{trl}(1:64,:) - mean(SL_raw.trial{corresponding_SL}(1:64,tw),2);
    end
    
    %% For both trial types...

    bad_channels = {};
    for e = 1:2 % 1 = SL, 2 = RL
        
        if e == 1
            data = SL_raw;
        elseif e == 2
            data = RL_raw;
        end
    
        % Downsample (mainly for the ICA)
        cfg = [];
        cfg.resamplefs = 200;
        cfg.detrend = 'no';
        cfg.demean = 'no';
        cfg.trials = 'all';
        cfg.sampleindex = 'no'; % whether to add a channel that has the original sample index or not
        data = ft_resampledata(cfg, data);
    
        disp([num2str(length(data.trial)) ' trials detected for ' subject])
        
        % Correct any channel swaps
        if ~isempty(swapped_channels{subject})
            for i = 1:size(swapped_channels{subject},1)
                from_idx = [];
                to_idx = [];
                channel_idx = [];
                disp(['Inserting data from ' swapped_channels{subject}{i,1} ' to ' swapped_channels{subject}{i,2} ' for ' schar])
                for chan = 1:length(data.label)
                    if strmatch(data.label{chan},swapped_channels{subject}{i,1})
                        from_idx = chan;
                    elseif strmatch(data.label{chan},swapped_channels{subject}{i,2})
                        to_idx = chan;
                    end
                end
                for trl = 1:length(data.trial)
                     data.trial{trl}(to_idx,:) = data.trial{trl}(from_idx,:);
                end
            end
        end

        % Find noisy channels
        cfg = [];
        cfg.viewmode = 'vertical';
        ft_databrowser(cfg, data);
        
        % check the percentage of time each channel exceeds 200 uV
        loud_channels = [];
        for trl = 1:length(data.trial)
            loud_channels(1:64,trl) = mean(data.trial{1}(1:64,:) > 200 | data.trial{1}(1:64,:) < -200,2);
        end
        figure; imagesc(loud_channels); set(gca,'YTick',1:64); set(gca,'YTickLabels',data.label(1:64));
        
        if e == 1
            bad_channels = input('Any bad channels to interpolate? Type {''Cz'',''Fp1''}, or 0 to skip: ');
            bad_channels = bad_channels';
        end
        
        if ~isempty(bad_channels) && iscell(bad_channels)
            cfg = [];
            cfg.badchannel = bad_channels;
            cfg.neighbours = neighbours;
            cfg.layout = layout;
            data = ft_channelrepair(cfg,data);
        end

        % Artifact correction via ICA (ignore bad channels)
        channel = {'eeg','-M1','-M2'};
        if ~isempty(bad_channels) && iscell(bad_channels)
            for chan = 1:length(bad_channels)
                channel{end+1} = ['-' bad_channels{chan}];
            end
        end
        channel = ft_channelselection(channel, data);
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
        
        if ~isempty(reject_comp)
            cfg = [];
            cfg.component = reject_comp; % to be removed component(s)
            data = ft_rejectcomponent(cfg, comp, data);
        end
        
        % Visual artifact rejection
        cfg = [];
        cfg.method = 'summary';
        cfg.channel = channel;
        cfg.keeptrial = 'yes';
        
        data.trialinfo = [1:length(data.trialinfo)]';
        clean = ft_rejectvisual(cfg,data);

        % List the rejected trials
        artifact = clean.cfg.artfctdef.summary.artifact;
        rejected = [];
        for i = 1:size(artifact,1)
            rejected(i,1) = find(clean.sampleinfo(:,1) == artifact(i,1));
        end
        keep_trials = setdiff(1:length(clean.trial),rejected);
        
        if e == 1
            idx = find(T.Subject == subjects(subj));
            T.SL_Artifact(idx(rejected)) = 1;
        elseif e == 2
            idx = find(T.Subject == subjects(subj) & ~isnan(T.RL_Trial));
            T.RL_Artifact(idx(rejected)) = 1;
        end
        
        % Select the non-artifactual trials
        cfg = [];
        cfg.trials = keep_trials;
        clean = ft_selectdata(cfg,clean);
        
        % Convert to scalp current density
        cfg = [];
        cfg.method = 'spline';
        cfg.layout = layout;
        cfg.neighbours = neighbours;
        cfg.trials = 'all';
        
        scd = ft_scalpcurrentdensity(cfg,clean);

        % Save data
        if e == 1
            tag = 'SL';
        elseif e == 2
            tag = 'RL';
        end
        
        disp(['Saving data for ' schar '...'])
        save(fullfile(dir_save,[schar '_' tag '_comp.mat']),'comp'); % ICA components
        save(fullfile(dir_save,[schar '_' tag '_data.mat']),'data'); % data with ICA correction
        save(fullfile(dir_save,[schar '_' tag '_clean.mat']),'clean'); % data with ICA correction & visual artifact rejection
        save(fullfile(dir_save,[schar '_' tag '_scd.mat']),'scd'); % scalp current density
        
        save(fullfile(dir_save,'T_table.mat'),'T'); % for group
        
    end
end

%% Plot group averages

trial_type = 'SL'; % 'SL' or 'RL'
data_type = 'scd'; % 'scd' or 'clean'

filelist = dir(fullfile(dir_save,['*' trial_type '_' data_type '.mat']));
load(fullfile(dir_save,'T_table.mat'));

% SL - 4 conditions
group = {};
for subj = 1:length(filelist)
    
    subject = strsplit(filelist(subj).name,'_');
    subject = sscanf(subject{1},'S%2d'); % get subject number
    
    load(fullfile(dir_save,filelist(subj).name));
    if strmatch(trial_type,'SL')
        idx = T(T.Subject == subject & ~isnan(T.SL_Trial),:);
    else idx = T(T.Subject == subject & ~isnan(T.RL_Trial),:);
    end
    
    if ~isempty(strmatch(data_type,'scd'))
        D = scd;
    elseif ~isempty(strmatch(data_type,'clean'))
        D = clean;
    end
    con_idx = idx.Condition(D.trialinfo);
    acc_idx = idx.Acc(D.trialinfo);
    rt_idx = idx.RT(D.trialinfo);
    
    for con = 1:4

        disp(['Averaging data for ' num2str(subject) ', condition ' num2str(con) '...'])
        
        if strmatch(trial_type,'SL')
            trial_idx = find(con_idx == con & acc_idx == 1);
        elseif strmatch(trial_type,'RL')
            trial_idx = find(con_idx == con & acc_idx == 1 & rt_idx > 1);
        end
        
        cfg = [];
        cfg.trials = trial_idx;
        data = ft_selectdata(cfg,D);
        
        cfg = [];
        cfg.channel = {'eeg','-M1','-M2'};
        cfg.keeptrials = 'no';
        group{subj,con} = ft_timelockanalysis(cfg,data);
        
    end
end

% Make grand average
cfg = [];
grand = {};
for con = 1:4
    grand{con} = ft_timelockgrandaverage(cfg,group{:,con});
end

% Look at overall topography
figure
cfg = [];
cfg.layout = layout;
ft_multiplotER(cfg,grand{:});

%% Cluster-based permutation

cfg = [];
ultragrand = ft_timelockgrandaverage(cfg,grand{:});

timewindows = [0 .5; .5 1; 1 1.5; 1.5 2; 2 2.5; 2.5 3];
% timewindows = [1.5 2.5];
figure
for t = 1:size(timewindows,1)
    subplot(3,2,t)
    cfg = [];
    cfg.layout = layout;
    cfg.xlim = timewindows(t,:);
%     cfg.zlim = [-30 30];
    cfg.gridscale = 300;
    cfg.markersymbol = '.';
    cfg.markersize = 20;
    cfg.colormap = viridis;
    ft_topoplotER(cfg,ultragrand);
    colorbar
end

figure
cfg = [];
cfg.layout = layout;
ft_multiplotER(cfg,ultragrand);

% Create contrasts
neut_pred = {};
fear_pred = {};
for subj = 1:size(group,1)
    
    cfg = [];
    cfg.operation = 'x1-x2';
    cfg.parameter = 'avg';
    neut_pred{subj} = ft_math(cfg,group{subj,2},group{subj,1});
    fear_pred{subj} = ft_math(cfg,group{subj,4},group{subj,3});
    
end
cfg = [];
grand_neut_pred = ft_timelockgrandaverage(cfg,neut_pred{:});
grand_fear_pred = ft_timelockgrandaverage(cfg,fear_pred{:});

figure
cfg = [];
cfg.layout = layout;
ft_multiplotER(cfg,grand_neut_pred,grand_fear_pred);

% cfg = [];
% channel = {'eeg','-M1','-M2'};
% channel = {'CPz','Pz'};
channel = {'PO8','P6','P8','P10','P4'};
% channel = {'P5','P7','P9'};
% channel = {'AFz','FPz','Fz','FCz','Cz'};
% channel = {'AF7','F3','FC1'};

latency = [0 1];

cfg = [];
cfg.layout = layout;
cfg.neighbours = neighbours;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = .05;
cfg.clusterstatistic = 'maxsum';
cfg.tail = 0;
cfg.clustertail = 0;
cfg.numrandomization = 1000;

nSub = size(group,1);
cfg.design(1,1:2*nSub)  = [ones(1,nSub) 2*ones(1,nSub)];
cfg.design(2,1:2*nSub)  = [1:nSub 1:nSub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.latency = latency;
cfg.channel = channel;

cfg.avgoverchan = 'yes';
cfg.avgovertime = 'no';

stat = ft_timelockstatistics(cfg,neut_pred{:},fear_pred{:});

% Show cluster probabilities
figure
imagesc(stat.prob < .05)
set(gca,'YTick',1:length(channel))
set(gca,'YTickLabels',stat.label)
colorbar

% cfg = [];
% cfg.layout = layout;
% cfg.neighbours = neighbours;
% ft_clusterplot(cfg,stat);

% Look at covariate
cfg = [];
cfg.layout = layout;
cfg.neighbours = neighbours;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesregrT';
cfg.correctm = 'cluster';
cfg.clusteralpha = .05;
cfg.clusterstatistic = 'maxsum';
cfg.tail = 0;
cfg.clustertail = 0;
cfg.numrandomization = 1000;

cfg.latency = latency;
cfg.channel = channel;
if length(cfg.channel) < length(grand_neut_pred.label)
    cfg.avgoverchan = 'yes';
end
cfg.avgovertime = 'no';

cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable

design = [];
these_subjects = unique(T.Subject);
for subj = 1:length(these_subjects)
    idx = find(T.Subject == these_subjects(subj) & T.Acc == 1 & T.RT > .5);
    this_rt = T.RT(idx);
    idx(zscore(this_rt) > 3 | zscore(this_rt) < -3,:) = [];
    this_rt = T.RT(idx);
    this_con = T.Condition(idx);
    design(subj,1) = mean(this_rt(this_con == 2))-mean(this_rt(this_con == 1));
    design(subj,2) = mean(this_rt(this_con == 4))-mean(this_rt(this_con == 3));
end
    
cfg.design = design(:,1);
stat_neut_cov = ft_timelockstatistics(cfg,neut_pred{:});

% Show cluster probabilities
figure
imagesc(stat_neut_cov.prob < .05)
set(gca,'YTick',1:length(channel))
set(gca,'YTickLabels',channel)
colorbar
title('Neutral Covariate')

cfg.design = design(:,2);
stat_fear_cov = ft_timelockstatistics(cfg,fear_pred{:});

% Show cluster probabilities
figure
imagesc(stat_fear_cov.prob < .05)
set(gca,'YTick',1:length(channel))
set(gca,'YTickLabels',channel)
colorbar
title('Fearful Covariate')

%% Custom plots
% Plot specific channel(s) for group with error bars

channel = {'PO8','P6','P8','P10','P4'};
% channel = {'AFz','FPz','Fz','FCz','Cz'};
% channel = {'CPz'};

x_time = grand{1,1}.time;

colours = [63 10 109;       % EN = purple
           17 118 216;      % UN = blue
           242 15 121;      % EF = pink
           242 173 15]/256; % UF = orange

contrast = 'expectation'; % 'all', or 'expectation'
       
figure
switch contrast
    case 'all'
        this_group = grand;
        this_individual = group;
        for con = 1:4

            channel_idx = [];
            for chan = 1:length(channel)
                for c = 1:length(this_group{1}.label)
                    if strmatch(this_group{1}.label{c},channel{chan})
                        channel_idx = [channel_idx,c];
                    end
                end
            end
            
            d = [];
            for subj = 1:size(this_individual,1)
                d(subj,:,:) = squeeze(this_individual{subj,con}.avg(channel_idx,:));
            end
            
            d = squeeze(mean(d,2));

            mu = mean(this_group{con}.avg(channel_idx,:),1);
            sem = std(d,1)/sqrt(size(d,1));
            upper = mu + sem;
            lower = mu - sem;

            x = [x_time fliplr(x_time)];
            patch = fill(x, [upper fliplr(lower)],colours(con,:),'HandleVisibility','off');
            set(patch,'edgecolor','none');
            set(patch,'FaceAlpha',.25);

            hold on
            if any(con == [1 3]) % expected = solid
                plot(x_time,mu,'LineWidth',2,'Color',colours(con,:));
            else % unexpected = dashed
                plot(x_time,mu,'LineWidth',2,'Color',colours(con,:),'LineStyle','--');
            end

        end
        xlim([x_time(1) x_time(end)])
        legend({'EN','UN','EF','UF'})
    case 'expectation'
        
        these_colours = [mean(colours(1:2,:),1);
                         mean(colours(3:4,:),1)];
        
         channel_idx = [];
         for chan = 1:length(channel)
             for c = 1:length(grand_neut_pred.label)
                 if strmatch(grand_neut_pred.label{c},channel{chan})
                     channel_idx = [channel_idx,c];
                 end
             end
         end
            
        d = [];
        for subj = 1:length(neut_pred)
            d(1,subj,:) = mean(neut_pred{subj}.avg(channel_idx,:),1);
            d(2,subj,:) = mean(fear_pred{subj}.avg(channel_idx,:),1);
        end
                     
        mu = [];
        mu(1,:) = mean(grand_neut_pred.avg(channel_idx,:),1);
        mu(2,:) = mean(grand_fear_pred.avg(channel_idx,:),1);
        
        for i = 1:2
            
            sem = std(squeeze(d(i,:,:)),1)/sqrt(size(squeeze(d(i,:,:)),1));
            upper = mu(i,:) + sem;
            lower = mu(i,:) - sem;

            x = [x_time fliplr(x_time)];
            patch = fill(x, [upper fliplr(lower)],these_colours(i,:),'HandleVisibility','off');
            set(patch,'edgecolor','none');
            set(patch,'FaceAlpha',.25);

            hold on
            plot(x_time,mu(i,:),'LineWidth',2,'Color',these_colours(i,:));
        end
        
        xlim([x_time(1) x_time(end)])
        plot([x_time(1) x_time(end)],[0 0],'k--')
        legend({'Neutral Surprise','Fearful Surprise'})
        
        stat_idx = stat.time(find(stat.prob < .05));
        scatter(stat_idx,ones(length(stat_idx),1)*-10,5,'filled','k'); hold on
        
        stat_idx = stat.time(find(stat_neut_cov.prob < .05));
        scatter(stat_idx,ones(length(stat_idx),1)*-11.5,5,'filled','MarkerFaceColor',these_colours(1,:)); hold on
        
        stat_idx = stat.time(find(stat_fear_cov.prob < .05));
        scatter(stat_idx,ones(length(stat_idx),1)*-13,5,'filled','MarkerFaceColor',these_colours(2,:)); hold on
        
end
title(strjoin(grand{1}.label(channel_idx),'  '))
set(gca,'TickLength',[0 0])

%% Draw covariate

% these_colours = [mean(colours(1:2,:),1);
%                  mean(colours(3:4,:),1)];
% 
% channel_idx = [];
% for chan = 1:length(channel)
%     for c = 1:length(grand_neut_pred.label)
%         if strmatch(grand_neut_pred.label{c},channel{chan})
%             channel_idx = [channel_idx,c];
%         end
%     end
% end
% 
% d = [];
% for subj = 1:length(neut_pred)
%     d(1,subj,:) = mean(neut_pred{subj}.avg(channel_idx,:),1);
%     d(2,subj,:) = mean(fear_pred{subj}.avg(channel_idx,:),1);
% end
% 
% mu = [];
% mu(1,:) = mean(d(1,:,stat_neut_cov.prob < .05),3);
% mu(2,:) = mean(d(1,:,stat_fear_cov.prob < .05),3);
% 
% figure
% for i = 1:2
%     subplot(1,2,i)
%     scatter(design(:,i)',mu(i,:),'filled','MarkerFaceColor',these_colours(i,:),'MarkerFaceAlpha',.5,'MarkerEdgeColor','k'); hold on
%     f = polyfit(design(:,i)',mu(i,:),1);
%     refline(f(1),f(2)); hold on
%     ax = gca;
%     plot([0 0],ax.YLim,'k'); hold on
%     plot(ax.XLim,[0 0],'k');
%     if i == 1
%         title([num2str(round(stat.time(find(stat_neut_cov.prob < .05,1,'first')),3)) ' to ' num2str(round(stat.time(find(stat_neut_cov.prob < .05,1,'last')),3)) ' ms'])
%     elseif i == 2
%         title([num2str(round(stat.time(find(stat_fear_cov.prob < .05,1,'first')),3)) ' to ' num2str(round(stat.time(find(stat_fear_cov.prob < .05,1,'last')),3)) ' ms'])
%     end
%     ylabel('More positive EEG activity for surprise')
%     xlabel('Faster for Surprise <---------> Faster for Expected')
% end

