%% Preprocess EEG Data in Fieldtrip

clear all
close all
clc

dir_ft = 'D:\group_garrido\Jessica\fieldtrip-20190419';
addpath(dir_ft)
ft_defaults;

dir_eeg = 'D:\group_garrido\Jessica\bcfs_eeg_reanalysis\data\Exp2\EEG';
dir_behav = 'D:\group_garrido\Jessica\bcfs_eeg_reanalysis\data\Exp2\Behavioural';

dir_save = 'G:\Private\Experiments\2_Experiment2_CFS\2019\eeg_FT_results';

subjects = setdiff(1:33,[3 23]);

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

downsample_avgs = 100; % 0 for no, X for new sampling rate
trial_type = 'SL'; % 'SL' or 'RL'
data_type = 'clean'; % 'scd' or 'clean'
restrict_standards = 1; % 0 for no, 1 for standard right BEFORE deviant, 2 for standard in middle of train, 3 for standard right AFTER deviant


filelist = dir(fullfile(dir_save,['*' trial_type '_' data_type '.mat']));
load(fullfile(dir_save,'T_table.mat'));

if restrict_standards > 0
   
    R1 = zeros(size(T,1),1);
    R2 = zeros(size(T,1),1);
    R3 = zeros(size(T,1),1);
    
    T.stan_before = R1;
    T.stan_middle = R2;
    T.stan_after = R3;
    
    for con = 1:2
        
        if con == 1
            % EN
            this_stan = 1;
            this_dev = 4;
            other_stan = 3;
        elseif con == 2
            % EF
            this_stan = 3;
            this_dev = 2;
            other_stan = 1;
        end
    
        % %% 1 = standard right BEFORE deviant
        r1 = find(T.Condition == this_dev) - 1;
        r1(T.Condition(r1) ~= this_stan,:) = [];
        R1(r1) = 1;

        % %% 3 = standard right AFTER deviant
        r3 = find(T.Condition == this_dev) + 1;
        r3(T.Condition(r3) ~= this_stan,:) = [];
        R3(r3) = 1;

        % %% 2 = standard in middle of train
        r2_idx = find(T.Condition == this_dev);
        r2 = [];
        for i = 1:length(r2_idx)
            
            if i == 1
                this_train = find(T.Condition == this_stan,1,'first'):r2_idx(i)-1;
            elseif i == length(r2_idx)
                this_train = r2_idx+1:find(T.Condition == this_stan,1,'last');
            else
                this_train = r2_idx(i-1)+1:r2_idx(i)-1;
                if any(T.Condition(this_train) == other_stan)
                    this_train = this_train(find(T.Condition(this_train) == other_stan,1,'last')+1:end);
                end
            end
            
            if mod(length(this_train)/2,1) == 0
                if rand(1) < .5
                    r2 = [r2; this_train(length(this_train)/2)];
                else 
                    r2 = [r2; this_train(length(this_train)/2 + 1)];
                end
            else
                r2 = [r2; this_train(ceil(length(this_train)/2))];
            end
            
        end
        R2(r2,1) = 1;
        
    end
    
    T.stan_before = R1;
    T.stan_middle = R2;
    T.stan_after = R3;
    
%     find(sum([T.stan_before, T.stan_middle T.stan_after],2) > 1);

end


% SL - 4 conditions
group = {};
trial_count = [];
for subj = 1:length(filelist)
    
    subject = strsplit(filelist(subj).name,'_');
    subject = sscanf(subject{1},'S%2d'); % get subject number
    
    disp(['Subject ' num2str(subject)])
    
    load(fullfile(dir_save,filelist(subj).name));
    
    if ~isempty(strmatch(data_type,'scd'))
        D = scd;
    elseif ~isempty(strmatch(data_type,'clean'))
        D = clean;
    end
    
    if length(D.label) < 64
%         error(' ')
    end
    
    idx = [];
    for trl = 1:length(D.trialinfo)
        if ~isempty(strmatch(trial_type,'SL'))
            idx(trl,:) = find(T.Subject == subject & T.SL_Trial == D.trialinfo(trl,1));
        elseif ~isempty(strmatch(trial_type,'RL'))
            idx(trl,:) = find(T.Subject == subject & T.RL_Trial == D.trialinfo(trl,1));            
        end
    end
    idx = T(idx,:);
    
    for con = 1:4

        disp(['Averaging data for ' num2str(subject) ', condition ' num2str(con) '...'])
        
        if (con == 1 || con == 3) && restrict_standards > 0
            if restrict_standards == 1
                rs_idx = idx.stan_before;
            elseif restrict_standards == 2
                rs_idx = idx.stan_middle;
            elseif restrict_standards == 3
                rs_idx = idx.stan_after;
            end
        else rs_idx = ones(size(idx,1),1);
        end
        
        if strmatch(trial_type,'SL')
            trial_idx = find(idx.Condition == con & idx.Acc == 1 & idx.SL_Artifact == 0 & rs_idx);
        elseif strmatch(trial_type,'RL')
            trial_idx = find(idx.Condition == con & idx.Acc == 1 & idx.RL_Artifact == 0 & rs_idx);
        end
%         trial_idx = find(isnan(idx.RT) & idx.SL_Artifact == 0);
%         trial_count(subj,1) = length(trial_idx);
        
        if ~isempty(trial_idx)
            cfg = [];
            cfg.trials = trial_idx;
            data = ft_selectdata(cfg,D);

            cfg = [];
            cfg.channel = {'eeg','-M1','-M2'};
            cfg.keeptrials = 'no';
            group{subj,con} = ft_timelockanalysis(cfg,data);
            
            if length(group{subj,con}.label) < 64
%                 error(' ' )
            end

            if downsample_avgs ~= 0
                cfg = [];
                cfg.resamplefs = downsample_avgs;
                cfg.detrend = 'no';
                cfg.demean = 'no';
                cfg.trials = 'all';
                cfg.sampleindex = 'no'; % whether to add a channel that has the original sample index or not
                group{subj,con} = ft_resampledata(cfg, group{subj,con});
            end
        end
    end
end
% save(fullfile(dir_save,['avresults_ds' num2str(downsample_avgs) '_' data_type '_rs' num2str(restrict_standards) '_nobc.mat']),'group')

% Re-apply baseline-correction
if ~isempty(strmatch(trial_type,'SL'))
    for subj = 1:size(group,1)
        for con = 1:size(group,2)
            cfg = [];
            cfg.demean = 'yes';
            cfg.baselinewindow = [-.1 0];
            group{subj,con} = ft_preprocessing(cfg,group{subj,con});
        end
    end
else
    tmp = load(fullfile(dir_save,['avresults_ds' num2str(downsample_avgs) '_SL_' data_type '_rs' num2str(restrict_standards) '_notbc.mat']));
    SL_group = tmp.group;
    i = 0;
    for subj = 1:size(group,1)
        for con = 1:size(group,2)
            if subj == 11 && con == 1
                i = i+1;
            end
            samp_idx = [1 find(abs(SL_group{subj+i,con}.time - 0) == min(abs(SL_group{subj+i,con}.time - 0)))];
            this_bc = mean(SL_group{subj+i,con}.avg(:,samp_idx(1):samp_idx(2)),2);
            group{subj,con}.avg = group{subj,con}.avg - this_bc;
        end
    end
end

% Make grand average
cfg = [];
cfg.keepindividual = 'yes';
grand = {};
for con = 1:4
    grand{con} = ft_timelockgrandaverage(cfg,group{:,con});
end

% % Look at overall topography
% figure
% cfg = [];
% cfg.layout = layout;
% ft_multiplotER(cfg,grand{:});

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cfg = [];
% ultragrand = ft_timelockgrandaverage(cfg,grand{:});
cfg = [];
ultragrand = {};
for con = 1:4
    ultragrand{con} = ft_timelockgrandaverage(cfg,group{:,con});
end
ultragrand = ft_timelockgrandaverage(cfg,ultragrand{:});

% timewindows = [0 .25; .25 .5; .5 .75; .75 1; 1 1.25; 1.25 1.5; 1.5 1.75; 1.75 2];
timewindows = [0.5 1];
% timewindows = [-1 -.5; -.2 0; 0 0];
figure
for t = 1:size(timewindows,1)
%     subplot(2,3,t)
    cfg = [];
    cfg.layout = layout;
    cfg.xlim = timewindows(t,:);
    cfg.zlim = [-3 3];
    cfg.gridscale = 300;
    cfg.markersymbol = '.';
    cfg.markersize = 20;
%     cfg.marker = 'labels';
%     cfg.colormap = viridis;
    cfg.style = 'straight';
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
neut = {};
fear = {};
expec = {};
unexpec = {};
interaction = {};
for subj = 1:size(group,1)
    
    cfg = [];
    cfg.operation = 'x1-x2';
    cfg.parameter = 'avg';
    neut_pred{subj} = ft_math(cfg,group{subj,2},group{subj,1});
    fear_pred{subj} = ft_math(cfg,group{subj,4},group{subj,3});
    interaction{subj} = ft_math(cfg,neut_pred{subj},fear_pred{subj});
%     neut_oddball{subj} = ft_math(cfg,group{subj,2},group{subj,3});
%     fear_oddball{subj} = ft_math(cfg,group{subj,4},group{subj,1});

    cfg = [];
    cfg.operation = '(x1+x2)/2';
    cfg.parameter = 'avg';
    neut{subj} = ft_math(cfg,group{subj,1},group{subj,2});
    fear{subj} = ft_math(cfg,group{subj,3},group{subj,4});
    expec{subj} = ft_math(cfg,group{subj,1},group{subj,3});
    unexpec{subj} = ft_math(cfg,group{subj,2},group{subj,4});

end
cfg = [];
cfg.keepindividual = 'yes';
grand_neut_pred = ft_timelockgrandaverage(cfg,neut_pred{:});
grand_fear_pred = ft_timelockgrandaverage(cfg,fear_pred{:});
grand_interaction = ft_timelockgrandaverage(cfg,interaction{:});
grand_neut_oddball = ft_timelockgrandaverage(cfg,neut_oddball{:});
grand_fear_oddball = ft_timelockgrandaverage(cfg,fear_oddball{:});


figure
cfg = [];
cfg.layout = layout;
ft_multiplotER(cfg,grand_neut_pred,grand_fear_pred);

save(fullfile(dir_save,['avresults_ds' num2str(downsample_avgs) '_' trial_type '_' data_type '_rs' num2str(restrict_standards)]),...
    'group','grand','ultragrand','neut_pred','fear_pred','grand_neut_pred','grand_fear_pred','neut','fear','expec','unexpec','interaction')
% end
%% Stats



% load(fullfile(dir_save,'avresults_ds100_SL_scd_rs1.mat'));

cfg = [];
% channel = {'eeg','-M1','-M2'};
% channel = {
%            'P7','P9','PO7','P8','P10','PO8',...
% %            'FCz','FC1','FC2','FC3','FC4',...
%            'Cz','C1','C2','C3','C4',...
%            'CPz','CP1','CP2','CP3','CP4',...
%            'Pz','P1','P2','P3','P4',...
% %            'FC5','FC6','C5','C6','CP5','CP6'
%             'Iz','Oz','O1','O2','POz'
%            };
channel = {'Iz','Oz','O1','O2','PO8','PO9','P7','P8','P9','P10','Pz','P1','P2','P3','P4',...
    'CPz','CP1','CP2','CP3','CP4','Cz','C1','C2','C3','C4'};

% if ~isempty(strmatch(trial_type,'SL'))
    latency = [0 1.5];
% elseif ~isempty(strmatch(trial_type,'RL'))
%     latency = [ultragrand.time(1) 0];
% end

cfg.layout = layout;
cfg.neighbours = neighbours;
cfg.parameter = 'avg';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
% cfg.correcttail = 'prob';
cfg.numrandomization = 1000;
cfg.minnbchan = 1;

nSub = size(group,1);
cfg.design = [];
cfg.design(1,1:2*nSub)  = [ones(1,nSub) 2*ones(1,nSub)];
cfg.design(2,1:2*nSub)  = [1:nSub 1:nSub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.latency = latency;
cfg.channel = channel;

cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';

% stat = ft_timelockstatistics(cfg,neut{:},fear{:});
% stat = ft_timelockstatistics(cfg,expec{:},unexpec{:});
stat = ft_timelockstatistics(cfg,neut_pred{:},fear_pred{:});
% stat = ft_timelockstatistics(cfg,group{:,1},group{:,2});
% stat = ft_timelockstatistics(cfg,group{:,3},group{:,4});
% stat = ft_timelockstatistics(cfg,neut_oddball{:},fear_oddball{:});

% Show cluster probabilities
figure
subplot(1,2,1)
imagesc(stat.prob)
set(gca,'YTick',1:length(stat.label))
set(gca,'YTickLabels',stat.label)
colorbar
subplot(1,2,2)
imagesc(stat.mask)
set(gca,'YTick',1:length(stat.label))
set(gca,'YTickLabels',stat.label)

sig_chan = stat.label(sum(stat.mask > 0,2) > 0)
stat.time(sum(stat.mask > 0,1) > 0)

% cfg = [];
% cfg.layout = layout;
% cfg.neighbours = neighbours;
% ft_clusterplot(cfg,stat);

%% Look at covariate
cfg = [];
cfg.layout = layout;
cfg.neighbours = neighbours;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesregrT';
cfg.correctm = 'cluster';
% cfg.alpha = .05;
% cfg.clusterstatistic = 'maxsum';
% cfg.tail = 0;
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;

cfg.latency = latency;
cfg.channel = channel;

cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';

cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable

if ~isempty(strmatch(trial_type,'SL'))
    eeg_subjects = [1:22 24:33];
else
    eeg_subjects = [1:2 4:22 24:33];
end

design = [];
for subj = 1:length(neut_pred)
    S = eeg_subjects(subj);
    
    design = [design; S];
    
    idx = find(T.Subject == S & T.Acc == 1 & T.RT > .5);
    this_rt = T.RT(idx);
    idx(zscore(this_rt) > 3 | zscore(this_rt) < -3,:) = [];
    this_rt = T.RT(idx);
    this_con = T.Condition(idx);
    design(subj,1) = mean(this_rt(this_con == 2))-mean(this_rt(this_con == 1));
    design(subj,2) = mean(this_rt(this_con == 4))-mean(this_rt(this_con == 3));
    design(subj,3) = mean(this_rt);
    
end
    
cfg.design = design(:,1);
% stat_neut_cov = ft_timelockstatistics(cfg,neut_pred{:});
stat = ft_timelockstatistics(cfg,group{:});

% Show cluster probabilities
figure
subplot(1,2,1)
imagesc(stat_neut_cov.prob)
set(gca,'YTick',1:length(stat_neut_cov.label))
set(gca,'YTickLabels',stat_neut_cov.label)
colorbar
title('Neutral Covariate')
subplot(1,2,2)
imagesc(stat_neut_cov.mask)
set(gca,'YTick',1:length(stat_neut_cov.label))
set(gca,'YTickLabels',stat_neut_cov.label)


cfg.design = design(:,2);
stat_fear_cov = ft_timelockstatistics(cfg,fear_pred{:});

% Show cluster probabilities
figure
subplot(1,2,1)
imagesc(stat_fear_cov.prob)
set(gca,'YTick',1:length(stat_fear_cov.label))
set(gca,'YTickLabels',stat_fear_cov.label)
colorbar
title('Fearful Covariate')
subplot(1,2,2)
imagesc(stat_fear_cov.mask)
set(gca,'YTick',1:length(stat_fear_cov.label))
set(gca,'YTIckLabels',stat_fear_cov.label)

% group all together



%% Custom plots
% Plot specific channel(s) for group with error bars

plot_individual_subjects = true;

% channel = stat.label(sum(stat.mask,2) > 0);

% channel = stat.label(sum(stat.mask > 0,2) > 0)';
% channel = {
%             'P7','P8','P9','PO7','PO8',...
%            'Pz','P1','P2','P3','P4',...
%             'CPz','CP1','CP2','CP3','CP4',...
%             'Cz','C1','C2','C3','C4',...
% %             'FCz','FC1','FC2','FC3','FC4',...
% % %             'FC5','FC6','C5','C6','CP5','Cp6'
%             };
% channel = {'P7','P8','P9','PO7','PO8'};
channel = {'Pz','P1','P2','P3','P4','CPz','CP1','CP2','CP3','CP4','Cz','C1','C1','C2','C3','C4'};
% channel = {'FCz','FC1','FC2','FC3','FC4','FC5','FC6','Cz','C1','C2','C3','C4','C5','C5','C6','CPz','CP1','CP2','CP3','CP4','CP5','CP6'};
% channel = {'CPz','CP1','CP2','CP3','CP4','Pz','P1','P2','P3','P4'};    

% if ~isempty(strmatch(trial_type,'SL'))
%     load('avresults_ds100_SL_clean_missedtrials.mat');
% end

x_time = grand{1,1}.time;

colours = [143 75 255;       % EN = purple
           51 125 255;      % UN = blue
           229 26 26;      % EF = pink
           255 147 0]/256; % UF = orange

contrast = 'all'; % 'all', or 'interaction', or 'oddball', or 'parameters'

figure

channel_idx = [];
for chan = 1:length(channel)
    for c = 1:length(grand{1}.label)
        if strmatch(grand{1}.label{c},channel{chan})
            channel_idx = [channel_idx,c];
        end
    end
end

switch contrast
    case 'all'
        
        if ~isempty(strmatch(trial_type,'SL')) && ~plot_individual_subjects
            mu = squeeze(mean(missed_grand.individual(:,channel_idx,:),2));
            mu = mean(mu,1);
            plot(x_time,mu,'Color',[.5 .5 .5],'LineWidth',2.5,'LineStyle',':'); hold on
        end
        
        if plot_individual_subjects
            S = 1:size(group,1);
        else S = 1;
        end
        
        for psubj = 1:length(S)
            
            if length(S) ~= 1
                subplot(4,8,psubj)
            end
            
            for con = 1:4

                if length(S) == 1
                    d = [];
                    for subj = 1:size(group,1)
                        d(subj,:,:) = squeeze(grand{con}.individual(subj,channel_idx,:));
                    end
                    d = squeeze(mean(d,2));

                    mu = mean(d,1);
                    sem = std(d,1)/sqrt(size(d,1));
                    upper = mu + sem;
                    lower = mu - sem;

                    x = [x_time fliplr(x_time)];
                    patch = fill(x, [upper fliplr(lower)],colours(con,:),'HandleVisibility','off');
                    set(patch,'edgecolor','none');
                    set(patch,'FaceAlpha',.25);
                    hold on
                else
                    mu = group{psubj,con}.avg(channel_idx,:);
                    mu = mean(mu,1);
                end

                if any(con == [1 3]) % expected = solid
                    plot(x_time,mu,'LineWidth',2,'Color',colours(con,:));
                else % unexpected = dashed
                    plot(x_time,mu,'LineWidth',2,'Color',colours(con,:),'LineStyle','--');
                end
                hold on

            end
        
        xlim([x_time(1) x_time(end)])
%         legend({'EN','UN','EF','UF'})
        
        end
        
    case 'interaction'
        
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
        d(1,:,:) = mean(grand_neut_pred.individual(:,channel_idx,:),2);
        d(2,:,:) = mean(grand_fear_pred.individual(:,channel_idx,:),2);

        d = d(:,:,:);
        
        mu = mean(d,2);
        
        for i = 1:2
            
            sem = std(squeeze(d(i,:,:)),1)/sqrt(size(d,2));
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
        
        stat_idx = stat.time(find(sum(stat.prob < .05,1) > 0));
        scatter(stat_idx,ones(length(stat_idx),1)*-.55,5,'filled','k'); hold on
%         
%         stat_idx = stat.time(find(sum(stat_neut_cov.prob < .05,1) > 0));
%         scatter(stat_idx,ones(length(stat_idx),1)*-11.5,5,'filled','MarkerFaceColor',these_colours(1,:)); hold on
%         
%         stat_idx = stat.time(find(sum(stat_fear_cov.prob < .05,1) > 0));
%         scatter(stat_idx,ones(length(stat_idx),1)*-13,5,'filled','MarkerFaceColor',these_colours(2,:)); hold on

    case 'oddball'

        these_colours = [mean(colours(2,:),1);
                         mean(colours(3,:),1)];
        
         channel_idx = [];
         for chan = 1:length(channel)
             for c = 1:length(grand_neut_oddball.label)
                 if strmatch(grand_neut_oddball.label{c},channel{chan})
                     channel_idx = [channel_idx,c];
                 end
             end
         end
            
        d = [];
        d(1,:,:) = mean(grand_neut_oddball.individual(:,channel_idx,:),2);
        d(2,:,:) = mean(grand_fear_oddball.individual(:,channel_idx,:),2);

        d = d(:,:,:);
        
        mu = mean(d,2);
        
        for i = 1:2
            
            sem = std(squeeze(d(i,:,:)),1)/sqrt(size(d,2));
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
        legend({'Neutral Oddball','Fearful Oddball'})
        
    case 'parameters'

        line_alpha = .5;
        sig_time = x_time(sum(stat.mask) > 0);
        
        for con = 1:4

            d = [];
            for subj = 1:size(group,1)
                d(subj,:,:) = squeeze(grand{con}.individual(subj,channel_idx,:));
            end
            d = squeeze(mean(d,2));
            mu = mean(d,1);

            if any(con == [1 3]) % expected = solid
                P = plot(x_time,mu,'LineWidth',3,'Color',colours(con,:));
            else % unexpected = dashed
                P = plot(x_time,mu,'LineWidth',3,'Color',colours(con,:),'LineStyle',':');
            end
            P.Color(4) = line_alpha;
            hold on
            
            if any(con == [1 3]) % expected = solid
                plot(sig_time,mu(1,sum(stat.mask) > 0),'LineWidth',3,'Color',colours(con,:));
            else % unexpected = dashed
                plot(sig_time,mu(1,sum(stat.mask) > 0),'LineWidth',3,'Color',colours(con,:),'LineStyle',':');
            end
            
            hold on

        end
        
        xlim([x_time(1) x_time(end)])
        legend({'EN','UN','EF','UF'})
        
end
title(strjoin(grand{1}.label(channel_idx),'  '))
set(gca,'TickLength',[0 0])
ax = gca;
if ~isempty(strmatch(trial_type,'SL'))
    plot(repmat(nanmean(T.RT(T.Acc == 1)),2,1),ax.YLim,'k--','LineWidth',2)
else
    plot([0 0],ax.YLim,'k--','Linewidth',2);
end

  
% Get specific data at time point (Range)
this_time = [0 0];
these_samples = [find(abs(x_time-this_time(1)) == min(abs(x_time-this_time(1)))),...
                 find(abs(x_time-this_time(2)) == min(abs(x_time-this_time(2))))];
savedata = [];             
for con = 1:length(grand)
    tmp = grand{con}.individual(:,channel_idx,:);
    tmp = mean(tmp,2);
    savedata(:,con) = mean(tmp(:,:,these_samples(1):these_samples(2)),3);
end
[h p ci stats] = ttest(savedata(:,1),savedata(:,2))
[h p ci stats] = ttest(savedata(:,3),savedata(:,4))

%% Draw covariate

x_time = neut_pred{1}.time;
sig_chan = find(sum(stat.mask,2) > 0);
sig_time = x_time(find(sum(stat.mask) > 0));

channel_idx = [];
for chan = 1:length(sig_chan)
    this_sig = stat.label(sig_chan(chan));
    for c = 1:length(grand{1}.label)
        if strmatch(grand{1}.label{c},this_sig{1})
            channel_idx = [channel_idx,c];
        end
    end
end

[sorted,sort_idx] = sort(cfg.design);

if length(sorted) > length(group)
    
    % match cfg.design to actual data
    con_idx = repmat(1:4,length(group),1);
    con_idx = con_idx(:);
    subj_idx = repmat(1:length(group),1,4)';
    
    sorted_full_idx = [subj_idx(sort_idx) con_idx(sort_idx)];
    
    % Split in half
    these_colours = [linspace(colours(1,1),colours(4,1),length(sorted)/2)',...
        linspace(colours(1,2),colours(4,2),length(sorted)/2)',...
        linspace(colours(1,3),colours(4,3),length(sorted)/2)'];
    
    figure
    
    subplot(1,2,1)
    d1 = cell(1,4);
    for i = 1:length(sorted)/2
        this_con = sorted_full_idx(i,2);
        d1{this_con} = [d1{this_con}; mean(group{sorted_full_idx(i,1),this_con}.avg(channel_idx,:),1)];
    end
    mu = [];
    for c = 1:4
        mu(c,:) = mean(d1{c});
        upper = mu(c,:) + std(mu(c,:))/sqrt(size(mu,1));
        lower = mu(c,:) - std(mu(c,:))/sqrt(size(mu,1));
        x = [x_time fliplr(x_time)];
        patch = fill(x, [upper fliplr(lower)],colours(c,:),'HandleVisibility','off');
        set(patch,'edgecolor','none');
        set(patch,'FaceAlpha',.25);

        hold on
    end  
    for c = 1:4 
        plot(x_time,mu(c,:),'LineWidth',2,'Color',colours(c,:));
    end
    set(gca,'TickLength',[0 0])
    
    subplot(1,2,2)
    d2 = cell(1,4);
    for i = (length(sorted)/2)+1:length(sorted)
        ii = i-(length(sorted)/2);
        this_con = sorted_full_idx(i,2);
        d2{this_con} = [d2{this_con}; mean(group{sorted_full_idx(i,1),this_con}.avg(channel_idx,:),1)];
    end
    mu = [];
    for c = 1:4
        mu(c,:) = mean(d2{c});
        upper = mu(c,:) + std(mu(c,:))/sqrt(size(mu,1));
        lower = mu(c,:) - std(mu(c,:))/sqrt(size(mu,1));
        x = [x_time fliplr(x_time)];
        patch = fill(x, [upper fliplr(lower)],colours(c,:),'HandleVisibility','off');
        set(patch,'edgecolor','none');
        set(patch,'FaceAlpha',.25);

        hold on
    end  
    for c = 1:4 
        plot(x_time,mu(c,:),'LineWidth',2,'Color',colours(c,:));
    end
    set(gca,'TickLength',[0 0])

end
    

% % Split into groups
% G = 2; % how many groups
% g = ceil(length(neut_pred)/G);
% 
% these_colours = [72 17 191; 112 200 255]/256;
% 
% if ~separate_plots
%     these_colours = [linspace(these_colours(1,1),these_colours(2,1),G)',...
%                      linspace(these_colours(1,2),these_colours(2,2),G)',...
%                      linspace(these_colours(1,3),these_colours(2,3),G)'];
%     figure
%     cc = 0;
%     if g == length(neut_pred)
%         g = 1;
%     end
%     for i = 1:g:length(neut_pred)
%         d = [];
%         cc = cc + 1;
%         for s = i:i+(g-1)
%             d(s-i+1,:) = mean(neut_pred{sort_idx(s,1)}.avg(sig_chan,:));
%         end
%         if g ~= 1
%             d = mean(d);
%         end
%         P = plot(x_time,d,'LineWidth',2,'Color',these_colours(cc,:)); hold on
%         P.Color(4) = 1;
%     end
%     xlim([x_time(1) x_time(end)])
%     plot([x_time(1) x_time(end)],[0 0],'k--')
%     
% else
%     
%     these_colours = colours;
%     
%     figure;
%     g_idx = 1:g:length(neut_pred);
%     for f = 1:G
%        
%         subplot(1,2,f);
%         d = [];
%         max_g = g_idx(f)+g-1;
%         if max_g > length(neut_pred)
%             max_g = length(neut_pred);
%         end
%         for s = g_idx(f):max_g
%             for c = 1:4
%                 d(s,c,:) = mean(group{sort_idx(s,1),c}.avg(sig_chan,:));
%             end
%         end
%         mu = squeeze(mean(d,1));
%         upper = mu + squeeze(std(d))/sqrt(size(d,1));
%         lower = mu - squeeze(std(d))/sqrt(size(d,1));
%         
%         for c = 1:2
%             
%             x = [x_time fliplr(x_time)];
%             patch = fill(x, [upper(c,:) fliplr(lower(c,:))],these_colours(c,:),'HandleVisibility','off'); hold on
%             set(patch,'edgecolor','none');
%             set(patch,'FaceAlpha',.25);
%         end
%         for c = 1:2
% %             if any(c == [1 3])
%                 plot(x_time,mu(c,:),'LineWidth',2,'Color',these_colours(c,:)); hold on
% %             else
% %                 plot(x_time,mu(c,:),'LineWidth',2,'Color',these_colours(c,:),'LineStyle','--'); hold on
% %             end
%         end
%         
%         xlim([x_time(1) x_time(end)])
%         plot([x_time(1) x_time(end)],[0 0],'k--')
%         ylim([-1 .75])
%         title(['UN - EN = ' num2str(round(mean(cfg.design(sort_idx(g_idx(f):max_g,1))),2))])
%         set(gca,'TickLength',[0 0])
%         
%     end
%     
% end
