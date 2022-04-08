clear all
clc

%% Directories

dir_data = 'D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\preprocessed';
dir_imgs = 'D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\images';

subjects = [1:2 4:22 24:33]; % 1:33;
N = length(subjects);

cmap = [77, 242, 255;
        98, 174, 255;
        255, 180, 29;
        255, 68, 142]/255;
    
conditionlabels     = {'EN','UN','EF','UF'};
standardType        = {'all','first', 'last'};
baselineCorrect     = [true, false];
removeBad           = true;
filterData          = false; % whether notch filters have been applied or not

if filterData
    filterTag = 'filtered';
else
    filterTag = 'unfiltered';
end

dir_data = fullfile(dir_data,filterTag);
dir_results = fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag);

% Load behavioural data
output = preprocess_behav(2,false);
T = output.trialdata;
clear output

%% Convert epoched data from SPM to Fieldtrip

smoothing = 0.2; % smoothing window, in seconds

dir_spm = 'D:\Toolboxes\spm12';
addpath(dir_spm);
spm('defaults','eeg')
        
% get neighbours
addpath('D:\Toolboxes\fieldtrip-20191119\template\layout')
addpath('D:\Toolboxes\fieldtrip-20191119\template\neighbours')

cfg = [];
cfg.method = 'template';
cfg.layout = 'biosemi64.lay';
neighbours = ft_prepare_neighbours(cfg);

for s = 1:N
   
    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end
    
    disp('====================================')
    disp(schar)
    disp('====================================')
    
    % Empty variables
    trialcount = nan(length(standardType),4);
    dir_results = fullfile(dir_data,'newepoched');
    
    % Load data
    if subject==26
        filename = ['c' schar '_r1_epoched_stimlocked.mat']; % merged
    else
        filename = [schar '_r1_epoched_stimlocked.mat'];
    end
    D = spm_eeg_load(fullfile(dir_results,filename));

    badchannels = D.chanlabels(D.badchannels);

    % Re-reference
    S           = struct;
    S.D         = D;
    S.refchan   = 'average';

    D = spm_eeg_reref_eeg(S);

    % Convert to Fieldtrip
    FT = ftraw(D);

    % Interpolate bad channels
    cfg = [];
    cfg.channel = {'EEG','-M1','-M2'};
    FT = ft_selectdata(cfg,FT);

    cfg = [];
    cfg.method = 'weighted';
    cfg.badchannel = badchannels';
    cfg.neighbours = neighbours;
    FT = ft_channelrepair(cfg,FT);

    % Smooth data
    for trl = 1:length(FT.trial)
        for chan = 1:length(FT.label)
            FT.trial{trl}(chan,:) = movmean(FT.trial{trl}(chan,:),smoothing*FT.fsample);
        end
    end
    
    % Get trial info
    [trialinfo,nTrls] = getTrialInfo(D);
    
    trialinfo.badTrial = zeros(nTrls,1);
    if removeBad
        trialinfo.badTrial(D.badtrials) = 1;
    end
    
    % Match with behavioural log
    thisT = T(T.Subject==subject,:);
    if subject==1
        thisT = thisT(78:end,:); % EEG recording started late - missed first 77 trials 
    elseif subject==15 || subject==16 || subject==18 || subject==19 || subject==20 || subject==33
        thisT = thisT(1:end-1,:); % last EEG trial missing
    end

    match_conditions = thisT.Condition==trialinfo.stimulusType;
    match_rt = thisT.RT - trialinfo.RT;
    trialinfo.RT(isnan(thisT.RT)~=isnan(trialinfo.RT)) = thisT.RT(isnan(thisT.RT)~=isnan(trialinfo.RT)) + nanmean(trialinfo.RT-thisT.RT); % if the response trigger is missing, use the behavioural log
    if ~all(match_conditions) || any(abs(match_rt)>0.25) || any(isnan(thisT.RT)~=isnan(trialinfo.RT))
        error(['Behaviour does not match with MEG for ' schar])
    end

    % Baseline correct
    cfg                 = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.05 0];
    bFT                 = ft_preprocessing(cfg,FT);

    % Save data types
    for b = 1:length(baselineCorrect)

        if baselineCorrect(b)
            thisD = bFT;
        else
            thisD = FT;
        end

%         % Convert to SCD
%         cfg             = [];
%         cfg.neighbours  = neighbours;
% 
%         cfg.method      = 'finite';
%         scd_finite      = ft_scalpcurrentdensity(cfg,thisD);
% 
%         cfg.method      = 'spline';
%         scd_spline      = ft_scalpcurrentdensity(cfg,thisD);
% 
%         cfg.method      = 'hjorth';
%         scd_hjorth      = ft_scalpcurrentdensity(cfg,thisD);
% 
%         % Save
%         save(fullfile(dir_results,[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),...
%                         'thisD','scd_finite','scd_spline','scd_hjorth','trialcount','badchannels','trialinfo','thisT','-v7.3'); 

        % Save
        save(fullfile(dir_results,[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),...
            'thisD','trialcount','badchannels','trialinfo','thisT','-v7.3');

    end

    % (clean up)
    delete(fullfile(dir_results,['M' filename]))
    delete(fullfile(dir_results,['M' filename(1:end-4) '.dat']))

end

%% Make plots

restoredefaultpath
addpath('D:\Toolboxes\fieldtrip-master')
ft_defaults;

st = 3; % index in standardType
b = 1;  % index in baselineCorrect
datatype = 'orig'; % 'orig' or 'scd_finite' or 'scd_spline' or 'scd_hjorth'

all_SL = cell(N,4);
all_missed = cell(N,1);
all_RL = cell(N,4);
all_SL_RT = cell(N,3);
all_RL_RT = cell(N,3);
all_behav = [];
for s = 1:N
    
    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end
    
    disp(['Reading in response-locked data for ' schar '...'])
    
    % Load data
    [SL,RL,T,idx] = select_data(fullfile(dir_data,'newepoched',[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),neighbours,standardType{st},datatype);
    
    all_behav = [all_behav; T(idx,:)];

    % Select trials for each plot type
    for c = 1:4
        cfg = [];
        cfg.channels = {'EEG','-M1','-M2'};
        cfg.trials = find(T.Condition==c & idx);
        all_SL{s,c} = ft_timelockanalysis(cfg,SL);
        all_RL{s,c} = ft_timelockanalysis(cfg,RL);
    end
    
    cfg = [];
    cfg.trials = find(isnan(T.RT));
    if ~isempty(cfg.trials)
        all_missed{s,1} = ft_timelockanalysis(cfg,SL);
    end

    q = quantile(T.RT(idx),4);
    for r = 1:3
        cfg = [];
        cfg.channels = {'EEG','-M1','-M2'};
        cfg.trials = find(T.RT > q(r) & T.RT < q(r+1) & idx);
        all_SL_RT{s,r} = ft_timelockanalysis(cfg,SL);
        all_RL_RT{s,r} = ft_timelockanalysis(cfg,RL);
    end
    
    clear SL
    clear RL
    clear T
    clear idx
    
end

save(['D:\bCFS_EEG_Reanalysis\results\plotdata_' datatype '.mat'],'all_SL','all_RL','all_SL_RT','all_RL_RT','all_missed','all_behav');

grand_SL = cell(1,5);
grand_RL = cell(1,5);
for c = 1:5
    if c==5
        grand_SL{c} = ft_timelockgrandaverage([],all_SL{:});
        grand_RL{c} = ft_timelockgrandaverage([],all_RL{:});
    else
        grand_SL{c} = ft_timelockgrandaverage([],all_SL{:,c});
        grand_RL{c} = ft_timelockgrandaverage([],all_RL{:,c});
    end
end

grand_SL_RT = cell(1,3);
grand_RL_RT = cell(1,3);
for c = 1:3
    grand_SL_RT{c} = ft_timelockgrandaverage([],all_SL_RT{:,c});
    grand_RL_RT{c} = ft_timelockgrandaverage([],all_RL_RT{:,c});
end

% do robust averaging for the missed trials (due to low trial count for some subjects producing noisy averages)
missed = [];
for s = 1:N
    missed(s,:,:) = all_missed{s}.avg;
end

grand_missed = grand_SL{1};
for chan = 1:size(missed,2)
    grand_missed.avg(chan,:) = spm_robust_average(squeeze(missed(:,chan,:)),1);
end

missed_tc = zeros(N,1);
for s = 1:N
    missed_tc(s,1) = length(all_missed{s}.cfg.trials);
end

% Plot multiplots for stimulus-locked and response-locked data
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.linewidth = 1.2;
cfg.linecolor = [0 232 255;
        51 110 255;
        255 186 48;
        255 0 0]/255;

figure
cfg.xlim = [-0.1 3];
ft_multiplotER(cfg,grand_SL{1:4});

figure
cfg.xlim = [-2 0.1];
ft_multiplotER(cfg,grand_RL{1:4});

cfg.linecolor = [243, 140, 255 ;
        178, 71, 255 ;
        66, 0, 213]/255;

figure
cfg.xlim = [-0.1 3];
ft_multiplotER(cfg,grand_SL_RT{1:3});
colormap(colours(256,'viridis'))

figure
cfg.xlim = [-2 0.1];
ft_multiplotER(cfg,grand_RL_RT{1:3});
colormap(colours(256,'viridis'))

% Pick timepoints of interest for topoplots
% --- SL
tp = [0.75 2];

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.style = 'straight';
cfg.markersymbol  = '.';
cfg.markersize = 6;
cfg.markercolor = [0 0 0];
cfg.comment = 'no';
cmap = colours(256,'viridis');

thiszlim = grand_SL{5}.avg(~ismember(grand_SL{5}.label,{'M1','M2','Nz','SNz','LH','RH','LV','UV','Status'}),:);
cfg.zlim = [min(thiszlim(:)) max(thiszlim(:))]*.75;

figure
for i = 1:length(tp)
    subplot(1,length(tp),i)
    cfg.xlim = repmat(tp(i),1,2);
    ft_topoplotER(cfg,grand_SL{5});
    colormap(cmap)
end

% --- RL
tp = [-.1];

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.layout = 'biosemi64.lay';
cfg.style = 'straight';
cfg.markersymbol  = '.';
cfg.markersize = 6;
cfg.markercolor = [0 0 0];
cfg.comment = 'no';
cmap = colours(256,'viridis');

figure
for i = 1:length(tp)
    cfg.xlim = repmat(tp(i),1,2);
    ft_topoplotER(cfg,grand_RL{5});
    colormap(cmap)
end

% Pick electrodes of interest and make pretty plots
thischan = cell(1,2);
thischan{1} = {'P7','P9','PO7','P8','P10','PO8'};
thischan{2} = {'C1','C2','P1','P2','CP1','CP2','Cz','CPz','Pz'};

cmap = [0 232 255;
        51 110 255;
        255 186 48;
        255 0 0;
        125 125 125]/255;

% --- SL
for chan = 1:length(thischan)
    
    figure
    for c = 1:5
        if c==5
            x = grand_missed.time;
            y = grand_missed.avg(ismember(grand_missed.label,thischan{chan}),:);
        else
            x = grand_SL{c}.time;
            y = grand_SL{c}.avg(ismember(grand_SL{c}.label,thischan{chan}),:);
        end
%         sem = grand_SL{c}.sem(ismember(grand_SL{c}.label,thischan{chan}),:);
        if size(y,1)>1 % more than one channel
            y = mean(y);
%             sem = mean(sem);
        end
%         patch([x fliplr(x)],[y+sem fliplr(y-sem)],cmap(c,:),'facealpha',.2,'edgealpha',0); hold on
        plot(x,y,'color',cmap(c,:),'linewidth',1.3); hold on
    end
    set(gca,'ticklength',[0 0])
    xlim(grand_SL{5}.time([1 end]))
    plot(grand_SL{5}.time([1 end]),[0 0],'k:','linewidth',1); hold on
    ax = gca;
    plot([0 0],ax.YLim,'k:','linewidth',1); hold on
    for c = 1:4
        plot(repmat(mean(all_behav.RT(all_behav.Condition==c)),2,1),ax.YLim,'color',cmap(c,:),'linestyle','--','linewidth',1.4); hold on
    end
end

% --- RL
thischan = {'Cz','CPz','Pz'};

figure
for c = 1:4
    x = grand_RL{c}.time;
    y = mean(grand_RL{c}.avg(ismember(grand_RL{c}.label,thischan),:));
%     sem = mean(grand_RL{c}.sem(ismember(grand_RL{c}.label,thischan),:));
%     patch([x fliplr(x)],[y+sem fliplr(y-sem)],cmap(c,:),'facealpha',.2,'edgealpha',0); hold on
    plot(x,y,'color',cmap(c,:),'linewidth',1.3); hold on
end
xlim([-1 .5])
set(gca,'ticklength',[0 0])
plot(grand_RL{5}.time([1 end]),[0 0],'k:','linewidth',1); hold on
ax = gca;
plot([0 0],ax.YLim,'k:','linewidth',1); hold on

%% Get correlation between channels and RT

st = 3; % standard type (1 = all, 2 = first, 3 = last)
b = 1; % baseline correction (1 = yes, 2 = no);
datatype = 'orig';

ST = cell(1,N);
for s = 1:N
   
    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end
    
    disp('====================================')
    disp(schar)
    disp('====================================')
            
    % Load
    [SL,~,T,idx] = select_data(fullfile(dir_data,'newepoched',[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),neighbours,standardType{st},datatype); 

    cfg = [];
    cfg.trials = find(idx);
    data = ft_selectdata(cfg,SL);

    % Get T maps
    cfg = [];
    cfg.neighbours = neighbours;
    cfg.channel = {'EEG','-M1','-M2'};

    cfg.method = 'montecarlo';
    cfg.tail = 0;
    cfg.alpha = .025;

    cfg.correctm = 'cluster';
    cfg.clusteralpha = .05;
    cfg.minnbchan = 2; % how many channels needed to form a cluster (minimum) - set to 0 to ignore this criterion
    cfg.numrandomization = 100;

    cfg.latency = [0 3];
    cfg.avgoverchan = 'no';
    cfg.avgovertime = 'no';

    cfg.statistic = 'ft_statfun_correlationT';
    cfg.ivar = 1;
    cfg.design = T.RT(idx)';

    stat = ft_timelockstatistics(cfg,data);

    % Save to group variable
    ST{s,1} = stat;
end


save(fullfile(dir_results,['rtcorr_' datatype '.mat']),'ST');

% Group level analysis
A = cell(N,1);
for s = 1:N
    
    thisstat = ST{s};
    A{s,1}.dimord = thisstat.dimord;
    A{s,1}.label = thisstat.label;
    A{s,1}.time = thisstat.time;
    A{s,1}.avg = zscore(thisstat.stat);

end

B = A;
for s = 1:N
    B{s}.avg = zeros(size(B{s}.avg,1),size(B{s}.avg,2));
end

cfg = [];
cfg.neighbours = neighbours;
cfg.channel = {'EEG','-M1','-M2'};

cfg.method = 'montecarlo';
cfg.tail = 0;
cfg.alpha = .025;

cfg.correctm = 'cluster';
cfg.clusteralpha = .05;
cfg.minnbchan = 0; % how many channels needed to form a cluster (minimum) - set to 0 to ignore this criterion
cfg.numrandomization = 500;

cfg.latency = [0 3];
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';

cfg.statistic = 'depsamplesT';
cfg.design = [ones(1,N) ones(1,N)*2; 1:N 1:N];
cfg.ivar = 1;
cfg.uvar = 2;

stat = ft_timelockstatistics(cfg,A{:},B{:});

sigtime = {};
sigchan = {};
cc = 0;
for posneg = 1:2

    if posneg==1 % positive clusters
        thiscluster = stat.posclusters;
        thismap = stat.posclusterslabelmat;
    elseif posneg==2 % negative clusters
        thiscluster = stat.negclusters;
        thismap = stat.negclusterslabelmat;
    end

    if ~isempty(thiscluster)
        nclusters = find(extractfield(thiscluster,'prob') <= .05);
        for n = 1:length(nclusters)
            cc = cc + 1;
            sigtime{cc} = stat.time(sum(thismap == n) > 0);
            sigchan{cc} = stat.label(sum(thismap == n,2) > 0);
        end
    end
end

for j = 1:length(sigtime)

    figure

    cfg = [];
    cfg.parameter = 'stat';
    cfg.layout = 'biosemi64.lay';

    cfg.xlim = sigtime{j}([1 end]);

    cfg.highlight = 'on';
    cfg.highlightchannel = sigchan{j};
    cfg.highlightsymbol = '.';
    cfg.highlightcolor = [1 1 1];
    cfg.highlightsize = 20;

    cfg.style = 'straight';
    cfg.markersymbol  = '.';
    cfg.markersize = 6;
    cfg.markercolor = [0 0 0];
    cfg.comment = 'no';

    ft_topoplotER(cfg,stat);
    colormap(colours(256,'viridis'))
    colorbar

    title([num2str(cfg.xlim(1)) ' to ' num2str(cfg.xlim(end))])
    caxis([-2.5 2.5])

end

% Extract significant channels per subject
sigchan = cell(N,1);
bestchan = cell(N,1);
cpp = [];
rt = [];
subidx = [];
for s = 1:N
    
    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end
    disp(['Getting CPP for subject ' schar '...'])
    
    thesechan = {};
    sigpos = find(extractfield(ST{s}.posclusters,'prob') < .05);
    signeg = find(extractfield(ST{s}.negclusters,'prob') < .05);
    
    if isempty(sigpos) %&& isempty(signeg)
        if ST{s}.posclusters(1).prob < ST{s}.negclusters(1).prob
            sigpos = 1;
        else
            signeg = 1;
        end
    end
    
    allbestchan = [];
    if ~isempty(sigpos)
        for i = 1:length(sigpos)
            
            sigmatrix = ST{s}.posclusterslabelmat==i;
            sigchanidx = sum(sigmatrix,2)>0;
            sigchanidx = find(sigchanidx & ismember(ST{s}.label,{'Cz','C1','C2','C3','C4',...
                                                          'Cpz','Cp1','Cp2','Cp3','Cp4',...
                                                          'Pz','P1','P2','P3','P4'}));
            sigchan{s} = ST{s}.label(sigchanidx);
            
            chanstats = nan(length(sigchanidx),1);
            for chan = 1:length(sigchanidx)
                 chanstats(chan,:) = mean(ST{s}.stat(chan,sigmatrix(sigchanidx(chan),:)));
            end
            [maxval,maxidx] = max(abs(chanstats));
            allbestchan = [allbestchan, sigchanidx(maxidx) maxval];
        end
    end
    if ~isempty(signeg)
        for i = 1:length(signeg)
            
            sigmatrix = ST{s}.negclusterslabelmat==i;
            sigchanidx = sum(sigmatrix,2)>0;
            sigchanidx = find(sigchanidx & ismember(ST{s}.label,{'Cz','C1','C2','C3','C4',...
                                                          'Cpz','Cp1','Cp2','Cp3','Cp4',...
                                                          'Pz','P1','P2','P3','P4'}));
            sigchan{s} = ST{s}.label(sigchanidx);
            
            chanstats = nan(length(sigchanidx),1);
            for chan = 1:length(sigchanidx)
                 chanstats(chan,:) = mean(ST{s}.stat(chan,sigmatrix(sigchanidx(chan),:)));
            end
            [maxval,maxidx] = max(abs(chanstats));
            allbestchan = [allbestchan, sigchanidx(maxidx) maxval];
        end
    end
    
    if ~isempty(allbestchan)
        bestchan{s} = ST{s}.label{allbestchan(allbestchan(:,2)==max(allbestchan(:,2)),1)};
    
        [SL,~,T,idx] = select_data(fullfile(dir_data,'newepoched',[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),neighbours,standardType{st},datatype); 
        fidx = find(idx);
        nTrls = sum(idx);

        thiscpp = nan(nTrls,length(SL.time{1}));
        for trl = 1:nTrls
             thiscpp(trl,:) = SL.trial{fidx(trl)}(ismember(SL.label,bestchan{s}),:);
        end

        [sorted sortidx] = sort(T.RT(idx));

        cpp = [cpp; thiscpp(sortidx,:)];
        rt = [rt; sorted];
        subidx = [subidx; repmat(s,size(thiscpp,1),1)];
    end
end
save(fullfile(dir_results,['rtcorr_' datatype '_cpp.mat']),'cpp','rt','subidx','sigchan','bestchan');

% Make all trials have the same polarity
peakamp = nan(size(cpp,1),1);
for trl = 1:size(cpp,1)
    peakamp(trl,1) = cpp(trl,round(rt(trl)*100)); 
end
cpp(peakamp<0,:) = cpp(peakamp<0,:) * (-1);

% Plot all
[sorted sortidx] = sort(rt);

figure
imagesc(imgaussfilt(normr(cpp(sortidx,:)),3)); hold on
colormap(colours(100,'viridis'));
caxis(caxis*.75)
plot(round(sorted*100),1:size(cpp,1),'k','linewidth',2)
set(gca,'ticklength',[0 0])

%% Do subject statistics

datatype = 'orig';
locking = 'response';
includeRT = false; % include trial-by-trial RT as a regressor at the first level
minclusternum = 2;

% load previous stats (if exists)
if includeRT
    savefilename = 'individualStats_RT';
else
    savefilename = 'individualStats';
end

savefilename = [savefilename '_' datatype];

switch locking
    case 'stimulus'
        savefilename = ['SL_' savefilename];
    case 'response'
        savefilename = ['RL_' savefilename];
end

savefilename = [savefilename '_minclusternum-' num2str(minclusternum)];

try
    load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,[savefilename '.mat']));
catch
    STATS = cell(1,N);
    interpchannels = cell(N,1);
end

for s = 1:N
   
    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end
    
    disp('====================================')
    disp(schar)
    disp('====================================')

    if ~includeRT
        ST = cell(length(standardType),length(baselineCorrect),9);
    else
        ST = cell(length(standardType),length(baselineCorrect),4);
    end
    if ~isempty(STATS{s})
        for st = 1:size(STATS{s},1)
            for b = 1:size(STATS{s},2)
                for i = 1:size(STATS{s},3)
                    if ~isempty(STATS{s}{st,b,i})
                        ST(st,b,i) = STATS{s}(st,b,i);
                    end
                end
            end
        end        
    end
    
    for st = 3%1:length(standardType)
        for b = 1%:length(baselineCorrect)
            
            % Load
            [SL,RL,T,idx] = select_data(fullfile(dir_data,'newepoched',[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),neighbours,standardType{st},datatype); 

            % If response-locked, make all trials the same length
            if strcmp(locking,'response')
                
                rtminmax = [-1 0];
  
                tlengths = nan(length(RL.trial),2);
                for trl = 1:length(RL.trial)
                    tlengths(trl,:) = RL.time{trl}([1 end]);
                end
                
                rmidx = tlengths(:,1) > rtminmax(1) | tlengths(:,2) < rtminmax(2);
                idx(rmidx) = 0;
                if any(rmidx & idx)
                    disp(['!!!!!!!!!! REMOVING ' num2str(sum(idx & rmidx)) ' TRIALS FOR BEING OUTSIDE OF RT LIMIT !!!!!!!!'])
                end
                
            end
            
            switch locking
                case 'stimulus'
                    data = cell(1,4);
                    cfg = [];
                    for c = 1:4
                        cfg.trials = find(idx & T.Condition==c);
                        data{c} = ft_selectdata(cfg,SL);
                    end
                case 'response'
                    data = cell(1,4);
                    cfg = [];
                    for c = 1:4
                        cfg.trials = find(idx & T.Condition==c);
                        data{c} = ft_selectdata(cfg,RL);
                    end
                    
                    cfg = [];
                    for c = 1:4
                        cfg.toilim = rtminmax;
                        data{c} = ft_redefinetrial(cfg,data{c});
                    end
            end

            %{
            cfg = [];
            cfg.layout = 'biosemi64.lay';
            cfg.linecolor = cmap;

            figure
            ft_multiplotER(cfg,data{:});
            %}
            
            % Get T maps
            cfg = [];
            cfg.neighbours = neighbours;
            cfg.channel = {'EEG','-M1','-M2'};

            cfg.method = 'montecarlo';
            cfg.tail = 0;
            cfg.alpha = .025;

            cfg.correctm = 'cluster';
            cfg.clusteralpha = .05;
            cfg.minnbchan = minclusternum;
            cfg.numrandomization = 100;

            switch locking
                case 'stimulus'
                    cfg.latency = [0 3];
                case 'response'
                    cfg.latency = [rtminmax(1) 0];
            end
            cfg.avgoverchan = 'no';
            cfg.avgovertime = 'no';

            % ------------------------

            if ~includeRT
                
                cfg.statistic = 'indepsamplesT';
                cfg.ivar = 1;

                for i = 1:12 %1:12 % emotion, expectation, (WITHIN BLOCK) interaction, neutral pred, fearful pred, (WITHIN EMOTION) interaction, neutral pred, fearful pred, EN, UN, EF, UF

                    if i==1 % neutral - fearful
                        
                        A = ft_appenddata([],data{1},data{2});
                        B = ft_appenddata([],data{3},data{4});
                        
                    elseif i==2 % expected - unexpected
                        
                        A = ft_appenddata([],data{1},data{3});
                        B = ft_appenddata([],data{2},data{4});
                        
                    elseif i==3 % within-block interaction effect

                        avEF = ft_timelockanalysis([],data{3});
                        avEN = ft_timelockanalysis([],data{1});
                        
                        A = data{2}; % UN
                        for trl = 1:length(A.trial)
                            A.trial{trl} = A.trial{trl} - avEF.avg; % minus average EF
                        end
                        
                        B = data{4}; % UF
                        for trl = 1:length(B.trial)
                            B.trial{trl} = B.trial{trl} - avEN.avg; % minus average EN
                        end

                    elseif i==4 % UN vs EF

                        A = data{2}; % UN
                        B = data{3}; % EF

                    elseif i==5 % UF vs EN

                        A = data{4}; % UF
                        B = data{1}; % EN
                        
                    elseif i==6 % within-emotion interaction effect

                        avEF = ft_timelockanalysis([],data{3});
                        avEN = ft_timelockanalysis([],data{1});
                        
                        A = data{2}; % UN
                        for trl = 1:length(A.trial)
                            A.trial{trl} = A.trial{trl} - avEN.avg; % minus average EN
                        end
                        
                        B = data{4}; % UF
                        for trl = 1:length(B.trial)
                            B.trial{trl} = B.trial{trl} - avEF.avg; % minus average EF
                        end

                    elseif i==7 % UN vs EN

                        A = data{2}; % UN
                        B = data{1}; % EN

                    elseif i==8 % UF vs EF

                        A = data{4}; % UF
                        B = data{3}; % EF

                    elseif i>=9 % all

                        A = ft_appenddata([],data{i-8});
                        
                        B = A;
                        for trl = 1:length(A.trial)
                            B.trial{trl} = zeros(size(B.trial{trl},1),size(B.trial{trl},2));
                        end

                    end

                    %{
                    pcfg = [];
                    pcfg.layout = 'biosemi64.lay';
                    pcfg.linecolor = cmap([1 4],:);
                    figure
                    ft_multiplotER(pcfg,ft_timelockanalysis([],A),ft_timelockanalysis([],B));
                    %}

                    cfg.design = [ones(1,length(A.trial)) ones(1,length(B.trial))*2];

                    stat = ft_timelockstatistics(cfg, A, B);
                    if any(~isreal(stat.stat(:)))
                        error(['Imaginary numbers detected for ' schar ', i = ' num2str(i)])
                    end
                    disp(['Lowest p value = ' num2str(min(stat.prob(:)))])

                    if any(isnan(stat.prob(:)))
                        error(['NaNs detected in stats for ' schar])
                    end

                    ST{st,b,i} = stat;

                end
            else
                
                avEF = ft_timelockanalysis([],data{3});
                avEN = ft_timelockanalysis([],data{1});

                A = data{2}; % UN
                for trl = 1:length(A.trial)
                    A.trial{trl} = A.trial{trl} - avEF.avg; % minus average EF
                end

                B = data{4}; % UF
                for trl = 1:length(B.trial)
                    B.trial{trl} = B.trial{trl} - avEN.avg; % minus average EN
                end
                
                rtA = thisT.RT(idx & thisT.Condition==2) - mean(thisT.RT(idx & thisT.Condition==3));
                rtB = thisT.RT(idx & thisT.Condition==4) - mean(thisT.RT(idx & thisT.Condition==2));
                
                for i = 1:3 % interaction, UN, UF

                    if i==1
                        cfg.statistic = 'depsamplesregrT';
                        cfg.ivar = 1;
                        cfg.uvar = 2;
                
                        cfg.design = [rtA' rtB'; ones(1,length(rtA')) ones(1,length(rtB'))*2];
                        stat = ft_timelockstatistics(cfg, A, B);
                    elseif i==2
                        cfg.statistic = 'indepsamplesregrT';
                        cfg.ivar = 1;
                        cfg = rmfield(cfg,'uvar');
                        
                        cfg.design = rtA';
                        stat = ft_timelockstatistics(cfg, A);
                    elseif i==3
                        cfg.statistic = 'indepsamplesregrT';
                        cfg.ivar = 1;
                        
                        cfg.design = rtB';
                        stat = ft_timelockstatistics(cfg, B);
                    end

                    disp(['Lowest p value = ' num2str(min(stat.prob(:)))])

                    if any(isnan(stat.prob(:)))
                        error(['NaNs detected in stats for ' schar])
                    end

                    ST{st,b,i} = stat;

                end
            end
        end
    end
    
    if any(~isreal(stat.stat(:)))
        error(['Imaginary numbers detected for ' schar ', i = ' num2str(i)])
    end

    STATS{s} = ST;
    
end

save(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,[savefilename '.mat']),...
    'STATS','-V7.3');

%% Group level

smoothing = []; % in ms ([] to skip smoothing)
datatype = 'orig';
locking = 'response'; % 'stimulus' or 'response'
subj_minclusternum = 2;
group_minclusternum = 2;
            
cfg = [];
cfg.method = 'template';
cfg.layout = 'biosemi64.lay';
neighbours = ft_prepare_neighbours(cfg);

if includeRT
    filename = ['individualStats_RT_' datatype];
else
    filename = ['individualStats_' datatype];
end

switch locking
    case 'stimulus'
        filename = ['SL_' filename];
    case 'response'
        filename = ['RL_' filename];
end

filename = [filename '_minclusternum-' num2str(subj_minclusternum)];
load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,[filename '.mat'])); % loads 'STATS' variable

load('D:\bCFS_EEG_Reanalysis\data\Exp2\behav\stai.mat'); % loads 'stai' variable
stai = table2array(struct2table(stai));
stai = stai(subjects,:);
for i = 1:size(stai,2)
    stai(:,i) = stai(:,i) - mean(stai(:,i));
end
stai = array2table(stai,'variablenames',{'full','state','trait'});

% drift diffusion parameters
ddm = readtable('D:\bCFS_EEG_Reanalysis\results\ddmparams_nozscore.csv'); 
ddm.Subject(ddm.Experiment==2) = ddm.Subject(ddm.Experiment==2)-31; % 31 subjects in experiment 1
ddm = ddm(ddm.Experiment==2,:);
ddm = ddm(ismember(ddm.Subject,subjects),:);

for st = 3%1:length(standardType)
    for b = 1%:length(baselineCorrect)
        stat = cell(1,9);
        for i = 1:9 % emotion, expectation, (WITHIN BLOCK) interaction, neutral pred, fearful pred, (WITHIN EMOTION) interaction, neutral pred, fearful pred, all conditions
            
            if strcmp(locking,'response')
                minrt = nan(N,3); % min time, max time, num samples
                for s = 1:N
                    minrt(s,:) = [STATS{s}{st,b,i}.time([1 end]) length(STATS{s}{st,b,i}.time)];
                end
                minrt = max(minrt(:,1));
            end
            
            % Plot the prediction error for each emotion type
            if i==3 || i==6 || i==9
                
                if i==9
                    C = 4;
                else
                    C = 2;
                end

                pdata = cell(N,2);
                for s = 1:N
                    for c = 1:C
                        if i<9
                            thisstat = STATS{s}{st,b,i+c};
                        else
                            thisstat = STATS{s}{st,b,i-1+c};
                        end
                        pdata{s,c}.dimord = thisstat.dimord;
                        pdata{s,c}.label = thisstat.label;
                        pdata{s,c}.time = thisstat.time;
                        
                        thisavg = thisstat.stat;
                        if i==6
                            thisavg = thisavg - thisavg(:,1);
                        end

%                         thisavg = zscore(thisavg')'; % z-score each channel across time
%                         for chan = 1:size(thisavg,1)
%                             thisavg(chan,:) = smooth(thisavg(chan,:),0.25*1024);
%                         end

                        pdata{s,c}.avg = thisavg;
                    end
                end
                
                gpdata = cell(1,C);
                for c = 1:C
                    gpdata{c} = ft_timelockgrandaverage([],pdata{:,c});
                end
                for c = 1:C % upper
                    gpdata{length(gpdata)+1} = gpdata{c};
                    gpdata{length(gpdata)}.avg = gpdata{length(gpdata)}.avg+(gpdata{length(gpdata)}.var/2);
                end
                for c = 1:C % lower
                    gpdata{length(gpdata)+1} = gpdata{c};
                    gpdata{length(gpdata)}.avg = gpdata{length(gpdata)}.avg-(gpdata{length(gpdata)}.var/2);
                end
                
                cmap = [77, 242, 255;
                    98, 174, 255;
                    255, 180, 29;
                    255, 68, 142]/255;


                figure
                cfg = [];
                cfg.layout = 'biosemi64.lay';
                if i==9
                    cfg.linecolor = cmap;
                else
                    cfg.linecolor = repmat(cmap([2 4],:),3,1);
                end
                cfg.linewidth = 1.4;
                cfg.linestyle = repmat({'-','--','--'},C,1);
                cfg.linestyle = cfg.linestyle(:);
                ft_multiplotER(cfg,gpdata{:});
                
            end
            
            % Get the t maps for this statistic
            if i<9
                A = cell(N,1);
            elseif i==9
                A = cell(N,4);
            end
            for s = 1:N

                if i<9
                    thisstat = STATS{s}{st,b,i};
                    A{s,1}.dimord = thisstat.dimord;
                    A{s,1}.label = thisstat.label;
                    A{s,1}.time = thisstat.time;

                    thisavg = thisstat.stat;
%                     if i==6
%                         thisavg = thisavg - thisavg(:,1);
%                     end
%                     thisavg = zscore(thisavg')'; % z-score each channel across time
%                     for chan = 1:size(thisavg,1)
%                         thisavg(chan,:) = smooth(thisavg(chan,:),0.25*1024);
%                     end

                    A{s,1}.avg = thisavg;
                elseif i==9
                    for c = 1:4
                        thisstat = STATS{s}{st,b,8+c};
                        A{s,c}.dimord = thisstat.dimord;
                        A{s,c}.label = thisstat.label;
                        A{s,c}.time = thisstat.time;
    
                        thisavg = thisstat.stat;
%                         thisavg = zscore(thisavg')'; % z-score each channel across time
%                         for chan = 1:size(thisavg,1)
%                             thisavg(chan,:) = smooth(thisavg(chan,:),0.25*1024);
%                         end
                        A{s,c}.avg = thisavg;
                    end

%                     % z-score across conditions
%                     for chan = 1:length(A{s,1}.label)
%                         for tp = 1:length(A{s,1}.time)
%                             
%                         end
%                     end

                end
                
                % For response-locking, make all subjects have the same length
                if strcmp(locking,'response')
                    if i==9
                        C = 4;
                    else
                        C=1;
                    end
                    for c = 1:C
                        thisx = findMin(minrt(1),A{s,c}.time):findMin(0,A{s,c}.time);
                        A{s,c}.time =  A{s,c}.time(thisx);
                        A{s,c}.avg = A{s,c}.avg(:,thisx);
                    end
                end
            end

            % Make dummy variable
            if i<9
                B = A;
                for s = 1:N
                    B{s}.avg = zeros(size(B{s}.avg,1),size(B{s}.avg,2));
                end
            end
            
            % Do t test
            cfg = [];
            cfg.neighbours = neighbours;
            cfg.channel = {'EEG','-M1','-M2'};

%             cfg.method = 'montecarlo';
%             cfg.tail = 0;
%             cfg.alpha = .025;
% 
%             cfg.correctm = 'cluster';
%             cfg.clusteralpha = .05;
%             cfg.minnbchan = group_minclusternum;
%             cfg.numrandomization = 500;

            cfg.method = 'analytic';
            cfg.correctm = 'fdr';
            cfg.tail = 0;
            cfg.alpha = 0.025;

            if i==9
                
%                 ddmvals = [ddm.boundary_C1 ddm.boundary_C2 ddm.boundary_C3 ddm.boundary_C4];
                ddmvals = [ddm.drift_C1 ddm.drift_C2 ddm.drift_C3 ddm.drift_C4];
%                 ddmvals = [ddm.nondecision_C1 ddm.nondecision_C2 ddm.nondecision_C3 ddm.nondecision_C4];
                
                x = cell(N,4); % eeg data
                y = nan(N,4); % ddm parameters
                for s = 1:size(A,1)
                    if any(ddm.Subject==subjects(s))
                        for c = 1:4
                            x{s,c} = A{s,c};
                            y(s,c) = ddmvals(ddm.Subject==subjects(s),c);
                        end
                    end
                end
                ridx = all(isnan(y),2);
                x = x(~ridx,:);
                y = y(~ridx,:);

                A = x(:);
                ddmvals = y(:);

                cfg.statistic = 'ft_statfun_correlationT'; % 'indepsamplesregrT'
                cfg.design = ddmvals - mean(ddmvals);
                cfg.ivar = 1;
            else
                cfg.statistic = 'indepsamplesT'; % 'indepsamplesregrT'
                cfg.design = [ones(1,N) ones(1,N)*2];
                cfg.ivar = 1;
%                 cfg.uvar = 2;
            end

            % ------------------------

            switch locking
                case 'stimulus'
                    cfg.latency = [0 3];
                case 'response'
                    cfg.latency = [minrt(1) 0];
            end

            cfg.avgoverchan = 'no';
            cfg.avgovertime = 'no';

            if i<9
                stat{i} = ft_timelockstatistics(cfg, A{:}, B{:});
            elseif i==9
                stat{i} = ft_timelockstatistics(cfg, A{:});
            end

            figure
            imagesc(stat{i}.mask); %imagesc(1-stat{i}.prob);
            colormap('hot')
            xlabel('Time')
            ylabel('Channels')
            set(gca,'ytick',1:length(stat{i}.label));
            set(gca,'yticklabels',stat{i}.label);
            title(['Lowest p = ' num2str(min(stat{i}.prob(:)))])
            set(gcf,'position',[440  -141 484 957])
            drawnow;
            
            % Show significant time window/sensors
            sigtime = {};
            sigchan = {};
            incluster = false;
            maskidx = sum(stat{i}.mask)>0;
            if any(maskidx)
                maskidx = [find(diff([0 maskidx])>0)', find(diff([maskidx 0])<0)'];
                for g = 1:size(maskidx,1) 
                    sigtime{g} = stat{i}.time(maskidx(g,1):maskidx(g,2));
                    sigchan{g} = stat{i}.label(sum(stat{i}.mask(:,maskidx(g,1):maskidx(g,2)),2)>0);
                end
            end

%             % Show significant clusters
%             sigtime = {};
%             sigchan = {};
%             cc = 0;
%             for posneg = 1:2
%                 
%                 if posneg==1 % positive clusters
%                     thiscluster = stat{i}.posclusters;
%                     thismap = stat{i}.posclusterslabelmat;
%                 elseif posneg==2 % negative clusters
%                     thiscluster = stat{i}.negclusters;
%                     thismap = stat{i}.negclusterslabelmat;
%                 end
%                 
%                 if ~isempty(thiscluster)
%                     nclusters = find(extractfield(thiscluster,'prob') <= .05);
%                     for n = 1:length(nclusters)
%                         cc = cc + 1;
%                         sigtime{cc} = stat{i}.time(sum(thismap == n) > 0);
%                         sigchan{cc} = stat{i}.label(sum(thismap == n,2) > 0);
%                     end
%                 end
%             end
%                
            for j = 1:length(sigtime)
                
                figure
                
                cfg = [];
                cfg.parameter = 'stat';
                cfg.layout = 'biosemi64.lay';

                cfg.xlim = sigtime{j}([1 end]);

                cfg.highlight = 'on';
                cfg.highlightchannel = sigchan{j};
                cfg.highlightsymbol = '.';
                cfg.highlightcolor = [1 1 1];
                cfg.highlightsize = 20;
                
                cfg.style = 'straight';
                cfg.markersymbol  = '.';
                cfg.markersize = 6;
                cfg.markercolor = [0 0 0];
                cfg.comment = 'no';

                ft_topoplotER(cfg,stat{i});
                colormap(colours(256,'viridis'))
                colorbar
                title([num2str(cfg.xlim(1)) ' to ' num2str(cfg.xlim(end))])
                
            end

            % plot t maps
            pdata = [];
            for s = 1:length(A)
                pdata(s,:,:) = A{s}.avg;
            end
            for j = 1:length(sigtime)
                
                figure

                if i==9 % sort by DDM parameter value

                    groups = quantile(ddmvals,3); % split into low/medium/high
                    y = [];
                    for c = 1:length(groups)-1
                        y(:,:,c) = squeeze(mean(pdata(ddmvals>=groups(c) & ddmvals<=groups(c+1),ismember(A{1}.label,sigchan{j}),:),2)); % average over electrodes
                    end
                    cmap = [102, 0, 173;
                        255, 162, 0]/255;

                else
                    y = squeeze(mean(pdata(:,ismember(A{1}.label,sigchan{j}),:),2)); % average over electrodes
                    cmap = [0 0 0];
                end

                for c = 1:size(y,3)
                    thisy = squeeze(y(:,:,c));
                    m = mean(thisy);
                    sem = std(thisy)/sqrt(size(thisy,1));
                    upper = m+sem;
                    lower = m-sem;
                    x = A{1}.time;
                    sigx = x(ismember(A{1}.time,sigtime{j}));
    
                    patch([x fliplr(x)],[upper fliplr(lower)],cmap(c,:),'facealpha',.2,'edgecolor','none'); hold on
                    plot(x,m,'color',cmap(c,:),'linewidth',1.4); hold on
                    plot(x([1 end]),[0 0],'color',cmap(c,:),'linestyle',':','linewidth',1.2); hold on
                    set(gca,'ticklength',[0 0])
                    if c==1
                        ax = gca;
                        scatter(sigx,repmat(ax.YLim(1),length(sigx),1),'markerfacecolor','k','markeredgecolor','none'); hold on
                    end
                end
            end
        end
    end
end

%% Fit curves to response-locked trials

% cmap = [98, 202, 255
%                 108, 43, 255 
%                 255, 196, 29 
%                 255, 29, 29 ]/255;
%             
% theseSubjects = setdiff(subjects,[3]); % subject 3 missed about half the trials
% thisN = length(theseSubjects);
% 
% analyse_channels = {'Cz','C1','C2';
%                     'CPz','CP1','CP2';
%                     'Pz','P1','P2'};
% 
% slopes = [];
% for s = 1:thisN
% 
%     subject = theseSubjects(s);
%     if subject < 10
%         schar = ['S0' num2str(subject)];
%     else
%         schar = ['S' num2str(subject)];
%     end
%     
%     % Get response-locked data
%     tmp = load(fullfile(dir_data,'newepoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc0.mat']));
%     thisD = tmp.thisD;
%     idx = tmp.idx;
%     thisT = tmp.thisT;
%     T = thisT(tmp.idx,:);
%     
%     if size(thisT,1) ~= length(thisD.trial)
%         error('Mismatch with behavioural file')
%     end
%     
%     % Check for bad channels
%     channels = find(~ismember(thisD.label,...
%         {'M1','M2','Nz','SNz','LH','RH','LV','UV','Status'}));
%     
%     nTrls = length(thisD.trial);
%     nChan = length(channels);
%     
%     d = nan(nChan,length(thisD.time{1}),nTrls);
%     for trl = 1:nTrls
%         d(:,:,trl) = thisD.trial{trl}(channels,:);
%     end
%     
%     d = std(reshape(d,size(d,1),size(d,2)*size(d,3)),[],2); % get standard deviation of each channel
%     
%     outliers = gesd(d,.005); % relatively more strict alpha criterion for what makes a channel an outlier
%     
%     badchannels = thisD.label(channels(outliers));
%     thisinterp{st,b} = badchannels;
%     
%     disp(['Interpolating: ' badchannels'])
%     
%     % Interpolate
%     cfg = [];
%     cfg.method = 'weighted';
%     cfg.badchannel = badchannels;
%     cfg.neighbours = neighbours;
%     interp = ft_channelrepair(cfg,thisD);
%     
%     % Shift time axis for response-locked data
%     nTrls = length(interp.trial);
%     offset = nan(nTrls,1);
%     for trl = 1:nTrls
%         offset(trl,1) = -findMin(thisT.RT(trl),interp.time{trl});
%     end
% 
%     cfg = [];
%     cfg.offset = offset;
%     interp = ft_redefinetrial(cfg,interp);
% 
%     % only select channels & trials of interest (i.e., those around CPz)
%     cfg = [];
%     cfg.channel = analyse_channels(:);
%     cfg.trials = find(tmp.idx);
%     SL = ft_selectdata(cfg,interp);
% 
%     nTrls = length(SL.trial);
%     
%     % collate data & average over channel groups & baseline-correct
%     Y = [];
%     X = [];
%     for trl = 1:nTrls
%         
%         y = SL.trial{trl};
%         x = SL.time{trl};
%         
% %         xidx = x >= x(1)+.1 & x <= 0.1;
%         xidx = x >= (-min(T.RT)-.1) & x <= 0.1;
%         
%         % baseline-correct
%         thisbc = mean(y(:,xidx & x >= 0),2);
%         y = y - thisbc;
%         
%         yidx = zeros(1,length(x));
%         yidx((length(yidx)-sum(xidx)+1):length(yidx)) = 1;
%         yidx = yidx==1;
%         
%         if size(analyse_channels,1) > 1
%             for i = 1:size(analyse_channels,1)
%                 y(i,:) = mean(y(ismember(SL.label,analyse_channels(i,:)),:)); 
%             end
%             y(i+1:end,:) = [];
%         end
%         
%         if isempty(X)
%             X = nan(1,length(x));
%             X(yidx) = x(xidx);
%         else
%             X(trl,:) = nan(1,length(x));
%             X(trl,yidx) = x(xidx);
%         end
%         
%         if isempty(Y)
%             Y = nan(1,size(y,1),size(y,2));
%             Y(1,:,yidx) = y(:,xidx);
%         else
%             Y(trl,:,yidx) = y(:,xidx);
%         end
%     end
%     
%     figure
%     set(gcf,'position',[139 74 1575 920])
%     sgtitle(schar)
%     cc = 0;
%     for chan = 1:length(analyse_channels)
%         for c = 1:4
%             
%             cc = cc + 1;
%             thisy = squeeze(Y(T.Condition==c,chan,:));
%             
%             origx = find(~isnan(max(X)));
%             
%             t0 = find(origx,1,'first');
%             t1 = find(origx,1,'last');
%             
%             thisx = max(X(:,origx));
%             avy = nanmean(thisy(:,origx));
% 
%             % smooth data using up to 10th-order polynomial
%             smoothavy = findpolyfit(thisx,avy,10);
%             
%             % get peaks in smoothed data
%             peaks = [1 find(ischange(smoothavy,'linear','threshold',10)) length(thisx)];
% 
%             % plot
%             subplot(length(analyse_channels),4,cc)
%             plot(max(X),nanmean(thisy),'k:'); hold on
%             plot(thisx,smoothavy,'k','linewidth',1.4); hold on
%             scatter(thisx(peaks),smoothavy(peaks),70,'marker','x','markeredgecolor','k','linewidth',1.5); hold on
%             
%             % get linear trends between each peak
%             theseslopes = [];
%             ii = 0;
%             for i = 1:length(peaks)-1
%                     
%                 % find an above-zero rate of change between the peaks as the starting point
%                 thist0 = find(thisx >= thisx(peaks(i)) & thisx <= thisx(peaks(i+1)),1,'first');
% 
%                 if ~isempty(thist0)
% 
%                     % fit linear model
%                     idx = thist0:peaks(i+1);
%                     [thisfit,stats] = polyfit(thisx(idx),avy(idx),1);
%                     fitvals = polyval(thisfit,thisx(idx));
%                     mse = mean((avy(idx)-fitvals).^2);
%                     rsquared = 1 - stats.normr/norm(avy(idx)) - mse;
% 
%                     % add to table
%                     ii = ii + 1;
%                     theseslopes = [theseslopes; array2table([theseSubjects(s) c ii,...
%                         thisx(peaks(i)) thisx(peaks(i+1)) thisx(idx(1)) thisx(idx(end)),...
%                         smoothavy(peaks(i)) smoothavy(peaks(i+1)),...
%                         thisfit(1) mse rsquared],...
%                         'variablenames',{'subject','condition','number',...
%                         'peak1','peak2','t0','t1',...
%                         'smoothpeakamp1','smoothpeakamp2',...
%                         'slope','mse','r2'})];
% 
%                     % plot
%                     plot(thisx(idx),fitvals,'color',cmap(c,:),'linewidth',1.4); hold on;
%                 end
%             end
%             slopes = [slopes; theseslopes];
%         end
%     end
% end

%% Select CPP component

for s = 1:N
    
    subjslopes = slopes(slopes.subject==subjects(s),:);
    
    % find channel with largest positive peak near response
    for c = 1:4
        tmp = subjslopes(subjslopes.condition==c,:);
        peaks = unique([tmp.peak1; tmp.peak2]);
    end
    
end



% % Get best model & channel per subject
% allmodels = {'linear','exp1','exp2'};
% groupparam = [];
% groupchannels = {};
% groupmodels = {};
% groupfit = [];
% for s = 1:thisN
%     
%     tmp = squeeze(modelfits(s,:,:));
%     [I,J] = find(tmp==min(tmp(:)));
%     
%     groupchannels{s} = analyse_channels{I};
%     groupmodels{s} = allmodels{J};
%     groupfit(s) = tmp(I,J);
%     groupparam(s,:,:) = squeeze(paramdiff(s,I,:,:));
%     
% end
% 
% d = nan(thisN,4);
% idx = contains(groupmodels,'exp2');
% d(idx,:) = squeeze(sum(groupparam(idx,:,[2 4]),3));
% d(~idx,:) = [];
% % d(~idx,:) = squeeze(groupparam(~idx,:,1));
% 
% % remove subjects with extreme values
% d(any(abs(zscore(d)) > 2,2),:) = [];
% 
% figure
% for c = 1:4
%     jitterx = (rand(size(d,1),1)-.5)*.25 + c; 
%     scatter(jitterx,d(:,c)); hold on
%     scatter(c,mean(d(:,c)),'markerfacecolor','k','markeredgealpha',0); hold on
%     plot([c c],mean(d(:,c)) + [std(d(:,c))/sqrt(size(d,1)) -std(d(:,c))/sqrt(size(d,1))],'k'); hold on
% end
% 
% [~,p1] = ttest(d(:,1),d(:,2))
% [~,p2] = ttest(d(:,3),d(:,4))
