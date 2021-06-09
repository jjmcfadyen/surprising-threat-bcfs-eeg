%% Plots

clear all
clc

%% Parameters

filterTag = 'unfiltered';

subjects = [1:22 24:33];
N = length(subjects);

schar = cell(1,N);
for s = 1:N
    subject = subjects(s);
    if subject < 10
        schar{s} = ['S0' num2str(subject)];
    else
        schar{s} = ['S' num2str(subject)];
    end
end

% Exclude participants from plots
eegparticipants = [1:22 24:33]; % available EEG data (S23 is missing - corrupted file, and also many missed responses)
excludeSubjects = [3 5 23]; % too many missed trials
theseSubjects = find(~ismember(eegparticipants,excludeSubjects));
thisN = length(theseSubjects);

% Start up FieldTrip
addpath('D:\Toolboxes\fieldtrip-20191119')
ft_defaults

cfg = [];
cfg.method = 'template';
cfg.layout = 'biosemi64.lay';
neighbours = ft_prepare_neighbours(cfg);

%% Get data and downsample for easier use / visualisation purposes only

tmp = load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'SL_individualStats.mat'));

badchannels = cell(1,N);
for s = 1:N
    badchannels{s} = {tmp.interpchannels{s}{3,1}};
end
clear tmp

% Read in epoched data and sort into conditions
SL = cell(N,4);
RL = cell(N,4);
missed = cell(N,1);
for s = 1:N
   
    disp('====================================')
    disp(schar{s})
    disp('====================================')
    
    % =====================
    % STIMULUS-LOCKED
    % =====================
    
    % Load epoched stimulus-locked data
    tmp = load(fullfile('D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\preprocessed',filterTag,'4_epoched',...
                [schar{s} '_FT_epoched_lastStandards_bc1.mat']));
       
    % Downsample
    cfg = [];
    cfg.resamplefs = 100; % in Hz
    tmp.thisD = ft_resampledata(cfg,tmp.thisD);
    
    % Interpolate stimulus-locked data (baseline corrected)
    cfg = [];
    cfg.method = 'weighted';
    cfg.badchannel = badchannels{s}{1};
    cfg.neighbours = neighbours;
    SL_interp = ft_channelrepair(cfg,tmp.thisD);
         
    % Select conditions
    for c = 1:4
        cfg = [];
        cfg.trials = find(tmp.idx & tmp.thisT.Condition==c);
        SL{s,c} = ft_selectdata(cfg,SL_interp);
        SL{s,c}.rt = tmp.thisT.RT(tmp.idx & tmp.thisT.Condition==c);
    end
    
    % Get missed trials
    cfg = [];
    cfg.trials = find(~tmp.trialinfo.badTrial & isnan(tmp.thisT.RT));
    missed{s,1} = ft_selectdata(cfg,SL_interp);
    
    % =====================
    % RESPONSE-LOCKED
    % =====================
    
    % For response-locked data, first load the non-baseline-corrected stimulus-locked data
    tmp = load(fullfile('D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\preprocessed',filterTag,'4_epoched',...
        [schar{s} '_FT_epoched_lastStandards_bc0.mat']));
    
    % Downsample
    cfg = [];
    cfg.resamplefs = 100; % in Hz
    tmp.thisD = ft_resampledata(cfg,tmp.thisD);
    
    % Use 100 ms post-response to baseline correct
    for trl = 1:length(tmp.thisD.trial)
        thisrt = tmp.thisT.RT(trl);
        thistime = tmp.thisD.time{trl};
        thisbc = findMin(thisrt,thistime):findMin(thisrt+.1,thistime);
        tmp.thisD.trial{trl} = tmp.thisD.trial{trl} - mean(tmp.thisD.trial{trl}(:,thisbc),2);
    end
    
    % Interpolate response-time-shifted data
    cfg = [];
    cfg.method = 'weighted';
    cfg.badchannel = badchannels{s}{1};
    cfg.neighbours = neighbours;
    RL_interp = ft_channelrepair(cfg,tmp.thisD);
    
    % Shift so that time 0 is response
    nTrls = length(RL_interp.trial);
    offset = nan(nTrls,1);
    for trl = 1:nTrls
        offset(trl,1) = -findMin(tmp.thisT.RT(trl),RL_interp.time{trl});
    end

    cfg = [];
    cfg.offset = offset;
    RL_interp = ft_redefinetrial(cfg,RL_interp);

    % Remove responses too close to the end of the trial
    minmax = nan(nTrls,2);
    for trl = 1:nTrls
        minmax(trl,:) = RL_interp.time{trl}([1 end]); 
    end

    rmidx = minmax(:,2) < .1;
    tmp.idx(rmidx,:) = 0;

    minmax = [max(minmax(tmp.idx,1)) min(minmax(tmp.idx,2))];
    
    % Select conditions
    for c = 1:4
        cfg = [];
        cfg.trials = find(tmp.idx & tmp.thisT.Condition==c);
        RL{s,c} = ft_selectdata(cfg,RL_interp);
        RL{s,c}.rt = tmp.thisT.RT(tmp.idx & tmp.thisT.Condition==c);
    end
    
    clear tmp
    clear RL_interp
    clear SL_interp
          
end

% Save
disp('Saving data...')
save(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_SL.mat'),'SL','-v7.3');
save(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_RL.mat'),'RL','-v7.3');
save(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_missed.mat'),'missed','-v7.3');

clear SL
clear RL
clear missed
clear thisD
clear tmp

%% SL ERPs - by condition

% Load 
load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_SL.mat'));

% Timelock
SL_condition = cell(N,4);
for s = 1:N
    for c = 1:4
        SL_condition{s,c} = ft_timelockanalysis([],SL{s,c});
    end
end

% Save
save(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,['SL_condition.mat']),'SL_condition');

% Get grand average
gSL_condition = cell(1,4);
for c = 1:6
    gSL_condition{c} = ft_timelockgrandaverage([],SL_condition{theseSubjects,c}); 
end

% Get overall average
aSL_condition = ft_timelockgrandaverage([],gSL_condition{1:4});

% Plot topography
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.xlim = [1.5 3];

figure
ft_topoplotER(cfg,aSL_condition);
colormap(colours(100,'viridis'))
colorbar

% Plot conditions
cmap = [77, 242, 255;
        98, 174, 255;
        255, 180, 29;
        255, 68, 142]/255;

figure
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.linecolor = cmap;
cfg.linewidth = 1.4;
ft_multiplotER(cfg,gSL_condition{:});

% Plot standards vs deviants
figure
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.linecolor = cmap([3 2],:);
cfg.linewidth = 1.4;
ft_multiplotER(cfg,gSL_condition{[3 2]});
sgtitle('Fearful Block')

figure
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.linecolor = cmap([1 4],:);
cfg.linewidth = 1.4;
ft_multiplotER(cfg,gSL_condition{[1 4]});
sgtitle('Neutral Block')
    
%% RL ERPs

% Load 
load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_RL.mat'));

% Timelock
RL_condition = cell(N,4);
for s = 1:N
    for c = 1:4
        disp(['Timelocking subject ' num2str(s) ', condition ' num2str(c) '...'])
        RL_condition{s,c} = ft_timelockanalysis([],RL{s,c});
    end
end

% Save
disp('Saving response-locked averages...')
save(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,['RL_condition.mat']),'RL_condition');

% Get epoch lengths
rtlength = nan(N,2);
for s = 1:N
    thisrt = [RL{s,1}.rt; RL{s,2}.rt; RL{s,3}.rt; RL{s,4}.rt];
    rtlength(s,:) = [min(thisrt) max(thisrt)];
end

% Get grand average
gRL_condition = cell(1,4);
for c = 1:4
    gRL_condition{c} = ft_timelockgrandaverage([],RL_condition{theseSubjects,c}); 
end

% Get overall average
aRL_condition = ft_timelockgrandaverage([],gRL_condition{:});

% Plot topography
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.xlim = [0 0];

figure
ft_topoplotER(cfg,aRL_condition);
colormap(colours(100,'viridis'))
colorbar

% Plot
cmap = [77, 242, 255;
        98, 174, 255;
        255, 180, 29;
        255, 68, 142]/255;

figure
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.linecolor = cmap;
ft_multiplotER(cfg,gRL_condition{:});

%% SL ERPs - by response time

load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_SL.mat'));
load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_RT.mat'));

SL_rt = cell(N,3); % fast, medium, slow
for s = 1:N
    
    disp(['Processing subject ' num2str(s) '...'])
    
    d = ft_appenddata([],SL{s,1},SL{s,2},SL{s,3},SL{s,4});
    thisrt = [RT{s,2}; RT{s,3}; RT{s,4}; RT{s,5}];
    
    rq = quantile(thisrt,2);
    rtbins = nan(length(thisrt),1);
    rtbins(thisrt < rq(1)) = 1;
    rtbins(thisrt >= rq(1) & thisrt <= rq(2)) = 2;
    rtbins(thisrt > rq(2)) = 3;
    
    for c = 1:3
        cfg = [];
        cfg.trials = find(rtbins==c);
        SL_rt{s,c} = ft_timelockanalysis(cfg,d);
    end
end

% Save
disp('Saving stimulus-locked averages binned by response time...')
save(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,['SL_rt.mat']),'SL_rt');

% Grand average
gSL_rt = cell(1,3);
for c = 1:3
    gSL_rt{c} = ft_timelockgrandaverage([],SL_rt{:,c});
end

% Plot
cmap = [244, 176, 255;
        178, 136, 255;
        131, 164, 255 ]/255; 
figure
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.linecolor = cmap;
cfg.linewidth = 1.4;
ft_multiplotER(cfg,gSL_rt{:});

%% Missed ERPs

% Load
load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_missed.mat')); % 'missed' variable

% Average per subject
tc = zeros(N,1);
SL_missed = cell(1,N);
for s = 1:N
    if length(missed{s}.trial) > 1
        SL_missed{s} = ft_timelockanalysis([],missed{s});
        tc(s,1) = length(missed{s}.trial);
    end
end
SL_missed = SL_missed(~cellfun('isempty',SL_missed));

% Grand average
gSL_missed = ft_timelockgrandaverage([],SL_missed{:}); 

% Add main condition
load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,['SL_condition.mat'])); % 'SL_condition'

% Get grand average
gSL_condition = cell(1,6);
for c = 1:6
    gSL_condition{c} = ft_timelockgrandaverage([],SL_condition{theseSubjects,c}); 
end

% Get overall average
aSL_condition = ft_timelockgrandaverage([],gSL_condition{1:4});

% Plot all
cmap = [90, 197, 255;
        126, 126, 126]/255;

figure
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.linecolor = cmap;
cfg.linewidth = 1.6;
ft_multiplotER(cfg,aSL_condition,gSL_missed);

%% Mismatch responses

load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_SL.mat'));
load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,'plot_RT.mat'));

% Compute neural mismatch response
mismatch = cell(N,2);
standards = cell(N,2);
deviants = cell(N,2);
for s = 1:N
    for c = 1:2 % neutral mismatch, fearful mismatch
        
        disp('-----------------------------------------------------------------------')
        disp(['Subject ' num2str(s) ' of ' num2str(N) ', condition ' num2str(c) '...'])
        disp('-----------------------------------------------------------------------')
        
        if c==1
            cons = [2 3]; % UN, EF
        elseif c==2
            cons = [4 1]; % UF, EN
        end
        
        % get average standard that accompanies this deviant
        standards{s,c} = ft_timelockanalysis([],SL{s,cons(end)});
        
        % get deviant response per trial
        deviants{s,c} = SL{s,cons(1)};
        
        % get mismatch response of each deviant vs average standard
        mismatch{s,c} = SL{s,cons(1)};
        for trl = 1:length(mismatch{s,c}.trial)
            mismatch{s,c}.trial{trl} = mismatch{s,c}.trial{trl} - standards{s,c}.avg;
        end
        
        % reduce RAM
        SL{s,cons(1)} = 0;
        SL{s,cons(end)} = 0;
        
    end
end

clear SL
save(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,['SL_mismatch.mat']),...
    'mismatch','deviants','standards','-v7.3');

% Get average across trials
avmismatch = cell(N,2);
for s = 1:N
    for c = 1:2
        avmismatch{s,c} = ft_timelockanalysis([],mismatch{s,c});
    end
end

gavmismatch = cell(1,2);
for c = 1:2
    gavmismatch{c} = ft_timelockgrandaverage([],avmismatch{theseSubjects,c});
end

% Plot
cmap = [98, 174, 255;
        255, 68, 142]/255;

figure
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.linewidth = 1.4;
cfg.linecolor = cmap;
ft_multiplotER(cfg,gavmismatch{:});

% Group by response time
showMismatch = false; % true = mismatch waveform, false = deviant and avStandard

nBins = 2; % 1 = deviants are fastest, 2 = deviants are slowest
mrt = cell(N,2);
bindeviant = cell(N,2,nBins);
for s = 1:N
    for c = 1:2
        
        if c==1
            cons = [2 3]; % UN, EF
        elseif c==2
            cons = [4 1]; % UF, EN
        end
        
        avstandard = mean(RT{s,cons(end)+1});
        thisrt = RT{s,cons(1)+1} - avstandard; % deviant minus standard
        
        rtbins = zeros(length(thisrt),1);
        rtbins(thisrt <= 0) = 1; % negative numbers = faster deviant response
        rtbins(thisrt > 0)  = 2; % positive numbers = slower deviant response
        
        mrt{s,c} = [thisrt rtbins];
        
        % select EEG trials
        for b = 1:nBins
            cfg = [];
            cfg.trials = find(rtbins==b);
            if showMismatch
                bindeviant{s,c,b} = ft_timelockanalysis(cfg,mismatch{s,c});
            else
                bindeviant{s,c,b} = ft_timelockanalysis(cfg,deviants{s,c});
            end
        end
        
    end
end

gavbinmismatch = cell(2,nBins);
for c = 1:2
    for b = 1:nBins
        gavbinmismatch{c,b} = ft_timelockgrandaverage([],bindeviant{:,c,b});
    end
end

if ~showMismatch
    gavstandards = cell(2,1);
    for c = 1:2
        gavstandards{c,1} = ft_timelockgrandaverage([],standards{:,c});
    end
end

% Plot  
for c = 1:2
    
    if c==1
        cmap = [46 205 255;
                28 123 191]/255;
    elseif c==2
        cmap = [255, 142, 23 
                255, 38, 67]/255;
    end

    figure
    cfg = [];
    cfg.layout = 'biosemi64.lay';
    
    cfg.linestyle = {'-','--'}; % solid = faster for deviants, dotted = slower for deviants
    if ~showMismatch
        cfg.linestyle{3} = '-'; % standards
    end
    
    cfg.linecolor = cmap; % lighter for solid lines, darker for dotted lines
    if ~showMismatch
        cfg.linecolor(3,:) = [.5 .5 .5]; % standards
    end
    
    cfg.linewidth = 1.4;
    
    pdata = gavbinmismatch(c,:);
    if ~showMismatch
        pdata{3} = gavstandards{c};
    end
    ft_multiplotER(cfg,pdata{:});
    if c==1
        sgtitle('Neutral Mismatch')
    elseif c==2
        sgtitle('Fearful Mismatch');
    end
    
end

