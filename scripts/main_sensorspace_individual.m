clear all
clc

includeRT = false;
locking = 'stimulus'; % 'stimulus' or 'response'

%% Directories

dir_data = 'D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\preprocessed';
dir_imgs = 'D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\images';

subjects = [1:22 24:33]; % 1:33;
N = length(subjects);

theseSubjects = setdiff(subjects,[3]); % subject 3 missed about half the trials
thisN = length(theseSubjects);

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

%% Convert epoched data from SPM to Fieldtrip

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
    dir_results = fullfile(dir_data,'4_epoched');
    if ~exist(dir_results)
        mkdir(dir_results)
    end
    
    % Load data
    if subject==26
        filename = ['c' schar '_r1_epoched_fulltrial.mat']; % merged
    else
        filename = [schar '_r1_epoched_fulltrial.mat'];
    end
    D = spm_eeg_load(fullfile(dir_data,'4_epoched',filename));

    badchannels = D.chanlabels(D.badchannels);
    
    % Get trial info
    [trialinfo,nTrls] = getTrialInfo(D);
    
    trialinfo.badTrial = zeros(nTrls,1);
    if removeBad
        trialinfo.badTrial(D.badtrials) = 1;
    end
    
    thisT = T(T.Subject==subject,:);
    if subject==1
        thisT = thisT(78:end,:); % EEG recording started late - missed first 77 trials 
    elseif subject==18
        thisT = thisT(1:end-1,:); % last EEG trial missing
    end

    % Get T-maps per subject
    cd(dir_results)
    for st = 1%1:length(standardType)
        
        % Only select trials not marked as bad
        idx = trialinfo.badTrial==0;

        % Only select trials that have a response
        idx = idx & ~isnan(thisT.RT); % right response triggers missing, so use behavioural log instead

        % Only select trials with correct responses
        idx = idx & thisT.Acc;
        
        % Remove trials where RTs are too fast (< 500 ms)
        idx = idx & thisT.RT > .5;

        % Remove trials where the RT was > 3 SDs from the mean RT
        rtoutliers = find(idx);
        rtoutliers(:,2) = abs(zscore(thisT.RT(idx))) > 3;
        idx(rtoutliers(rtoutliers(:,2)==1)) = 0;
        
        % Select the type of standard
        switch standardType{st}
            case 'first'
                idx = idx & (contains(trialinfo.standardType,'deviant') | contains(trialinfo.standardType,'first'));
            case 'last'
                idx = idx & (contains(trialinfo.standardType,'deviant') | contains(trialinfo.standardType,'last'));
        end

        for c = 1:4
            trialcount(st,c) = sum(idx & trialinfo.stimulusType==c);
        end
    
        for b = 1:length(baselineCorrect)

            if baselineCorrect(b)
                if ~exist(fullfile(dir_data,'4_epoched',['b' filename]))
                    S           = [];
                    S.D         = D;
                    S.timewin   = [D.time(1) 0]*1000; % in ms
                    thisD       = spm_eeg_bc(S);
                else
                    thisD = spm_eeg_load(fullfile(dir_data,'4_epoched',['b' filename]));
                end
            else
                thisD = D;
            end

            % Convert to fieldtrip
            thisD = ftraw(thisD);
            
            % Save
            save(fullfile(dir_results,[schar '_FT_epoched_' standardType{st} 'Standards_bc' num2str(baselineCorrect(b)) '.mat']),...
                'thisD','trialcount','badchannels','idx','trialinfo','thisT','-v7.3');
            
        end
    end
end

%% Make plots

all_SL = cell(thisN,4);
all_missed = cell(thisN,1);
all_RL = cell(thisN,4);
all_behav = [];
for s = 1:thisN
    
    subject = theseSubjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end
    
    disp(['Reading in response-locked data for ' schar '...'])
    
    % Load data & interpolate
    [SL,RL,T,idx] = clean_data(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc1.mat']));
    
    % Select behavioural data
    behav = array2table([T.Subject(idx) T.Condition(idx) T.RT(idx) nan(sum(idx),1)],'variablenames',{'Subject','Condition','RT','StandardType'});
    behav.StandardType = T.Type(idx);
    
    all_behav = [all_behav; behav];
    
    clear behav
    
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
    
    clear SL
    clear RL
    clear T
    clear idx
    
end

save('D:\bCFS_EEG_Reanalysis\results\plotdata.mat','all_SL','all_RL','all_missed','all_behav');

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

% do robust averaging for the missed trials (due to low trial count for some subjects producing noisy averages)
missed = [];
for s = 1:thisN
    missed(s,:,:) = all_missed{s}.avg;
end

grand_missed = grand_SL{1};
for chan = 1:size(missed,2)
    grand_missed.avg(chan,:) = spm_robust_average(squeeze(missed(:,chan,:)),1);
end

missed_tc = zeros(thisN,1);
for s = 1:thisN
    missed_tc(s,1) = length(all_missed{s}.cfg.trials);
end

% Plot multiplots for stimulus-locked and response-locked data
cfg = [];
cfg.layout = 'biosemi64.lay';

figure
ft_multiplotER(cfg,grand_SL{5});

figure
ft_multiplotER(cfg,grand_RL{5});

% Pick timepoints of interest for topoplots
% --- SL
tp = [0.75 2.25];

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.style = 'straight';
cfg.markersymbol  = '.';
cfg.markersize = 6;
cfg.markercolor = [0 0 0];
cfg.comment = 'no';
cmap = colours(100,'viridis');

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
tp = [1 2];

cfg = [];
cfg.layout = 'biosemi64.lay';
cmap = colours(100,'viridis');

figure
for i = 1:length(tp)
    cfg.xlim = repmat(tp(i),2,1);
    ft_topoplotER(cfg,all_SL{:});
    colormap(cmap)
end

% Pick electrodes of interest and make pretty plots
thischan = cell(1,2);
thischan{1} = {'P7','P9','PO7','P8','P10','PO8'};
thischan{2} = {'Cz','C1','C2','CPz','CP1','CP2','Pz','P1','P2'};

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
            y = mean(grand_missed.avg(ismember(grand_missed.label,thischan),:));
        else
            x = grand_SL{c}.time;
            y = mean(grand_SL{c}.avg(ismember(grand_SL{c}.label,thischan),:));
        end
    %     sem = mean(grand_SL{c}.sem(ismember(grand_SL{c}.label,thischan),:));
    %     patch([x fliplr(x)],[y+sem fliplr(y-sem)],cmap(c,:),'facealpha',.2,'edgealpha',0); hold on
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
thischan = {'Cz','C1','C2','CPz','CP1','CP2','Pz','P1','P2'};

figure
for c = 1:4
    x = grand_RL{c}.time;
    y = mean(grand_RL{c}.avg(ismember(grand_RL{c}.label,thischan),:));
%     sem = mean(grand_SL{c}.sem(ismember(grand_SL{c}.label,thischan),:));
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
    [SL,~,T,idx] = clean_data(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc' num2str(baselineCorrect(b)) '.mat'])); 

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
    cfg.minnbchan = 0; % how many channels needed to form a cluster (minimum) - set to 0 to ignore this criterion
    cfg.numrandomization = 100;

    cfg.latency = [0 3];
    cfg.avgoverchan = 'no';
    cfg.avgovertime = 'no';

    cfg.statistic = 'indepsamplesregrT';
    cfg.ivar = 1;
    cfg.design = T.RT(idx);

    stat = ft_timelockstatistics(cfg,data);

    % Save to group variable
    ST{s,1} = stat;
end



save(fullfile(dir_results,'rtcorr_subjectLevel.mat'),'ST');

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

figure
for j = 1:length(sigtime)

    subplot(1,length(sigtime),j)

    cfg = [];
    cfg.parameter = 'stat';
    cfg.layout = 'biosemi64.lay';

    cfg.xlim = sigtime{j}([1 end]);

    cfg.highlight = 'on';
    cfg.highlightchannel = sigchan{j};
    cfg.highlightsymbol = '.';
    cfg.highlightcolor = [1 1 1];

    cfg.style = 'straight';
    cfg.markersymbol  = '.';
    cfg.markersize = 6;
    cfg.markercolor = [0 0 0];
    cfg.comment = 'no';

    ft_topoplotER(cfg,stat);
    colormap(colours(100,'viridis'))
    colorbar

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
    
        [SL,~,T,idx] = clean_data(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc' num2str(baselineCorrect(b)) '.mat'])); 
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
save(fullfile(dir_results,'rtcorr_subjectLevel_cpp.mat'),'cpp','rt','subidx','sigchan','bestchan');

% Make all trials have the same polarity
peakamp = nan(size(cpp,1),1);
for trl = 1:size(cpp,1)
    peakamp(trl,1) = cpp(trl,round(rt(trl)*1024)); 
end
cpp(peakamp<0,:) = cpp(peakamp<0,:) * (-1);

% Plot all
[sorted sortidx] = sort(rt);

figure
imagesc(imgaussfilt(normr(cpp(sortidx,:)),5)); hold on
colormap(colours(100,'viridis'));
caxis(caxis/2)
plot(round(sorted*1024),1:size(cpp,1),'k','linewidth',2)
set(gca,'ticklength',[0 0])

% Plot CPP data sorted by RT
[sorted sortidx] = sort(rtamp.RT);

figure
imagesc(imgaussfilt(normr(cpp(sortidx,:)),2)); hold on
colormap(colours(100,'viridis'))
plot(round(sorted*1024),1:size(rtamp,1),'k','linewidth',3); hold on
xlabel('Time (samples')
ylabel('Trial')
title(['Average amplitude of ' channels ', sorted by RT'])


%% Do subject statistics

% restoredefaultpath
addpath('D:\Toolboxes\fieldtrip-20191119')
ft_defaults

cfg = [];
cfg.method = 'template';
cfg.layout = 'biosemi64.lay';
neighbours = ft_prepare_neighbours(cfg);

STATS = cell(1,N);
interpchannels = cell(N,1);
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
        ST = cell(length(standardType),length(baselineCorrect),6);
    else
        ST = cell(length(standardType),length(baselineCorrect),4);
    end
    thisinterp = cell(length(standardType),length(baselineCorrect));
    for st = 3%1:length(standardType)
        for b = 1%:length(baselineCorrect)
            
            % Load
            [SL,RL,T,idx] = clean_data(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc' num2str(baselineCorrect(b)) '.mat'])); 

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
            cfg.minnbchan = 0; % how many channels needed to form a cluster (minimum) - set to 0 to ignore this criterion
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

                for i = 1:3 % interaction, neutral, fearful

                    if i==1 % within-block

                        avEF = ft_timelockanalysis([],data{3});
                        avEN = ft_timelockanalysis([],data{1});
                        
                        A = data{2}; % UN
                        for trl = 1:length(A.trial)
                            A.trial{trl} = A.trial{trl} - avEN.avg; % minus average EF
                        end
                        
                        B = data{4}; % UF
                        for trl = 1:length(B.trial)
                            B.trial{trl} = B.trial{trl} - avEF.avg; % minus average EN
                        end

                    elseif i==2 % UN vs EF

                        A = data{2}; % UN
                        B = data{1}; % EF

                    elseif i==3 % UF vs EN

                        A = data{4}; % UF
                        B = data{3}; % EN

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
    
    STATS{s} = ST;
    interpchannels{s} = thisinterp;
    
end

if includeRT
    savefilename = 'individualStats_RT.mat';
else
    savefilename = 'individualStats.mat';
end

switch locking
    case 'stimulus'
        savefilename = ['SL_' savefilename];
    case 'response'
        savefilename = ['RL_' savefilename];
end

save(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,savefilename),'STATS','interpchannels');

%% SVD

% load(fullfile(dir_results,'rtcorr_subjectLevel_cpp.mat')); % 'cpp','rt','subidx','sigchan','bestchan'

gAmp = [];
gLat = [];
gDat = [];
gBF = [];
gVF = [];
gOrig = [];
allBetas = cell(1,N);
figure
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

    for st = 3%1:length(standardType)
        for b = 1%:length(baselineCorrect)
            
            % Load
            [SL,~,T,idx] = clean_data(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc' num2str(baselineCorrect(b)) '.mat'])); 
            
            fidx = find(idx);
            nTrls = length(fidx);

            % Filter at 6 Hz
            cfg = [];
            cfg.trials = fidx;
%             if ~isempty(sigchan{s})
%                 cfg.channel = sigchan{s};
%             else
%                 cfg.channel = SL.label(ismember(SL.label,{'Cz','C1','C2','CPz','CP1','CP2','Pz','P1','P2'}));
                cfg.channel = SL.label(~ismember(SL.label,{'M1','M2','Nz','SNz','LH','RH','LV','UV','Status'}));
%             end
            nChan = length(cfg.channel);
            cfg.lpfilter = 'yes';
            cfg.lpfreq = 6;
            SL = ft_preprocessing(cfg,SL);
            
            % Get average amplitude
            y = nan(nTrls,nChan,length(SL.time{1}));
            for trl = 1:nTrls
                 y(trl,:,:) = SL.trial{trl};
            end
            
            % SVD
            m = squeeze(mean(y)); % average over trials to get mean response overall
            [bf,vf] = svdyeahyouknowme(m',3,false);
            
            % (correct sign)
            bfscaled = bf(:,1) + bf(:,2)*[-1,1];
            [C,lags] = xcorr(bfscaled(:,1)',bfscaled(:,2)','coeff');
            lagshift = lags(abs(C)==max(abs(C)));
            if lagshift < 0 % if higher BF2 values shift the waveform to the left (i.e., earlier), then flip sign
                bf(:,2) = bf(:,2)*(-1);
            end
            
%             subplot(8,4,s)
%             P = plot(bf(:,1) + bf(:,2)*linspace(-1,1,10),'linewidth',1.4); hold on
%             cmap = colours(length(P),'viridis');
%             for i = 1:length(P)
%                 P(i).Color = cmap(i,:);
%             end
%             plot(bf(:,1),'k','linewidth',2.2)
%             title(schar);
%             drawnow
            
            gBF(s,:,:) = bf;
            gVF(s,:) = vf;
            
            % Apply weights to individual trial data
            newy = nan(nTrls,size(y,3));
            betas = nan(nTrls,size(gBF,3)); % amplitude, latency
            for trl = 1:nTrls
                tmp = squeeze(y(trl,:,:))';
                tmp = tmp - mean(tmp);
                B = pinv(bf)*tmp;
                yhat = bf*B;
                newy(trl,:) = yhat*vf;
                betas(trl,:) = mean(B,2);
            end
            allBetas{s} = T(idx,:);
            allBetas{s}.beta1 = betas(:,1);
            allBetas{s}.beta2 = betas(:,2);
            allBetas{s}.beta3 = betas(:,3);
            
            % save to group variable (average across trials)
            conditions = T.Condition(idx);
            for c = 1:4
                gAmp(s,c) = mean(betas(conditions==c,1));
                gLat(s,c) = mean(betas(conditions==c,2));
                gDat(s,c,:) = mean(newy(conditions==c,:));
                gOrig(s,c,:,:) = mean(y(conditions==c,:,:));
            end            
        end
    end
end
save(fullfile(dir_results,'svd_bf1-3.mat'),'gAmp','gLat','gDat','gBF','gVF','gOrig','allBetas');

% Basis function
y = squeeze(mean(gBF));
sem = squeeze(std(gBF))/sqrt(size(gBF,1));
upper = [y+sem]';
lower = [y-sem]';

x = linspace(-.1,3,size(gBF,2));
cmap = [0, 255, 97;
    188, 113, 255;
    255, 170, 0]/255;

figure
for i = 1:3
    patch([x fliplr(x)],[upper(i,:) fliplr(lower(i,:))],cmap(i,:),'facealpha',.2,'edgecolor','none','handlevisibility','off'); hold on
    plot(x,y(:,i),'color',cmap(i,:),'linewidth',1.4); hold on
end
xlim(x([1 end]))
set(gca,'ticklength',[0 0])
xlabel('Time (seconds)')
title('Temporal basis function')

% First basis predicting second basis
bf = squeeze(mean(gBF));

figure
for i = 1:2
    subplot(1,2,i)
    x = linspace(-0.1,3,length(bf));
    P = plot(x,bf(:,1) + bf(:,i+1)*linspace(-1,1,10),'linewidth',1.4); hold on
    cmap = colours(length(P),'viridis');
    for i = 1:length(P)
        P(i).Color = cmap(i,:);
    end
    plot(x,bf(:,1),'k','linewidth',2)  
    title(['First basis scaled by basis #' num2str(i+1)])
    xlim(x([1 end]))
    set(gca,'ticklength',[0 0])
end

% Betas
width = 0.25;
linewidth = 2.5;
cmap = [0, 242, 242;
    0, 135, 255;
    255, 181, 23;
    255, 0, 0]/255;

for i = 1:3 % basis functions
    
    % remove outlier trials per participant
    dat = [];
    for s = 1:N
        if i==1
            y = allBetas{s}.beta1;
        elseif i==2
            y = allBetas{s}.beta2;
        elseif i==3
            y = allBetas{s}.beta3;
        end
        outliers = abs(zscore(y)) > 3;
        disp(['Removing ' num2str(sum(outliers)) ' for subject ' num2str(s)])
        for c = 1:4
            dat(s,c) = mean(zscore(y(allBetas{s}.Condition==c & ~outliers)));
        end
    end
    
    figure
    for c = 1:4 % condition

        % points
        [x,y] = beeswarm(dat(:,c),4,width);
        scatter(x+c,y,'markeredgecolor',cmap(c,:),'markerfacecolor',cmap(c,:),'markerfacealpha',.25,'markeredgealpha',.5); hold on

        % boxplot
        patch([c-width c+width c+width c-width],[quantile(y,.75),quantile(y,.75),quantile(y,.25),quantile(y,.25)],cmap(c,:),...
            'facealpha',.5,'edgecolor',cmap(c,:)); hold on
        plot([c-width c+width],repmat(quantile(y,.5),2,1),'color',cmap(c,:),'linewidth',linewidth); hold on
        plot([c-width c+width],repmat(quantile(y,.25),2,1),'color',cmap(c,:),'linewidth',linewidth); hold on
        plot([c-width c+width],repmat(quantile(y,.75),2,1),'color',cmap(c,:),'linewidth',linewidth); hold on
        plot([c-width c-width],quantile(y,[.25 .75]),'color',cmap(c,:),'linewidth',linewidth); hold on
        plot([c+width c+width],quantile(y,[.25 .75]),'color',cmap(c,:),'linewidth',linewidth); hold on

    end
end

% LME
ldat = [];
for s = 1:N
    tmp = allBetas{s};
    tmp.beta1 = zscore(tmp.beta1);
    tmp.beta2 = zscore(tmp.beta2);
    tmp.beta3 = zscore(tmp.beta3);
    outliers = abs(tmp.beta1) > 3 | abs(tmp.beta2) > 3 | abs(tmp.beta3) > 3;
    ldat = [ldat; tmp(~outliers,:)];
end

lme = fitglme(ldat,'beta1 ~ Emotion*Expectation + (1|Subject)')
lme = fitglme(ldat,'beta2 ~ Emotion*Expectation + (1|Subject)')
lme = fitglme(ldat,'beta3 ~ Emotion*Expectation + (1|Subject)')


figure
for i = 1:3
    
    if i==1
        tmp = ldat.beta1;
    elseif i==2
        tmp = ldat.beta2;
    elseif i==3
        tmp = ldat.beta3;
    end
    
    y = [];
    for s = 1:N
        for c = 1:4
            y(s,c) = mean(tmp(ldat.Subject==subjects(s) & ldat.Condition==c));
        end
    end
    
    m = mean(y);
    sem = std(y)/sqrt(N);
    upper = m+sem;
    lower = m-sem;
    
    subplot(1,3,i)
    bar(m); hold on
    for c = 1:4
        plot([c c],[lower(c) upper(c)],'k');
    end
end


% Cross-correlation between conditions
conditiontype = 'all'; % 'all' (4 conditions) or 'mismatch-within' (EN-UN/EF-UF) or 'mismatch-between' (EN-UF/EF-UN)

allpred = [];
allvpred = [];
for s = 1:N
    
    % get predicted data
    thisbf = squeeze(gBF(s,:,:));
    pred = [];
    switch conditiontype
        case 'all'
            for c = 1:4
                dat = squeeze(gOrig(s,c,:,:))';
                B = pinv([thisbf ones(length(thisbf),1)]) * dat;
                pred(c,:,:) = thisbf * B(1:end-1,:);
            end
        case 'mismatch-within'
            for c = 1:2
                if c==1
                    dat = squeeze(gOrig(s,2,:,:))' - squeeze(gOrig(s,1,:,:))'; % UN - EN
                elseif c==2
                    dat = squeeze(gOrig(s,4,:,:))' - squeeze(gOrig(s,3,:,:))'; % UF - EF
                end
                B = pinv([thisbf ones(length(thisbf),1)]) * dat;
                pred(c,:,:) = thisbf * B(1:end-1,:);
            end
        case 'mismatch-between'
            for c = 1:2
                if c==1
                    dat = squeeze(gOrig(s,2,:,:))' - squeeze(gOrig(s,3,:,:))'; % UN - EF
                elseif c==2
                    dat = squeeze(gOrig(s,4,:,:))' - squeeze(gOrig(s,1,:,:))'; % UF - EN
                end
                B = pinv([thisbf ones(length(thisbf),1)]) * dat;
                pred(c,:,:) = thisbf * B(1:end-1,:);
            end
    end
    
    vpred = [];
    for c = 1:size(pred,1)
        vpred(c,:) = squeeze(pred(c,:,:)) * gVF(s,:)';
    end
    
    allpred(s,:,:,:) = pred;
    allvpred(s,:,:) = vpred;
end

x = linspace(-.1,3,size(allpred,3));
if strcmp(conditiontype,'all')
    cmap = [0, 242, 242;
        0, 135, 255;
        255, 181, 23;
        255, 0, 0]/255;
else
    cmap = [ 0, 135, 255;
    255, 0, 0]/255;
end

figure
for chan = 1:size(allpred,4)
    subplot(8,8,chan);
    title(SL.label{chan})
    for c = 1:size(allpred,2)
        y = squeeze(allpred(:,c,:,chan));
        m = mean(y);
        sem = std(y)/sqrt(N);
        upper = m+sem;
        lower = m-sem;
        
        patch([x fliplr(x)],[upper fliplr(lower)],cmap(c,:),'facealpha',.2,'edgecolor','none'); hold on
        plot(x,m,'color',cmap(c,:),'linewidth',1.4); hold on
    end
end

switch conditiontype
    case 'all'
        combos = [
            1 2 3 4 
            1 3 2 4
            1 1 2 2;
            1 1 3 3;
            3 3 4 4;
            2 2 4 4
            ];
    case {'mismatch-within','mismatch-between'}
        combos = [1 2];
end

startpoint = findMin(0.25,x); % only look from 250 ms onwards

xcoef = [];
% peaktime = []; % time of first peak
for s = 1:N

    pred = squeeze(allpred(s,:,:,:));
    
    % do cross-correlation for each channel
    thisxcoef = nan(size(combos,1),size(pred,3));
    for c = 1:size(combos,1)
        for chan = 1:size(pred,3)
            
            if strcmp(conditiontype,'all')
                y1 = mean(pred(combos(c,1:2),:,chan),1);
                y2 = mean(pred(combos(c,3:4),:,chan),1);
            else
                y1 = squeeze(pred(combos(c,1),:,chan));
                y2 = squeeze(pred(combos(c,2),:,chan));
            end
            [C,lags] = xcorr(y1,y2,'coef');
            thisxcoef(c,chan) = lags(abs(C)==max(abs(C)));
        end
    end
    xcoef(s,:,:) = thisxcoef;
    
%     for chan = 1:size(pred,3)
%         switch conditiontype
%             case 'all'
%                 y = squeeze(pred(:,:,chan));
%             case 'mismatch-within'
%                 y = [squeeze(pred(2,:,chan))-squeeze(pred(1,:,chan));   % UN-EN
%                      squeeze(pred(4,:,chan))-squeeze(pred(3,:,chan))];  % UF-EF
%             case 'mismatch-between'
%                 y = [squeeze(pred(2,:,chan))-squeeze(pred(3,:,chan));   % UN-EF
%                      squeeze(pred(4,:,chan))-squeeze(pred(1,:,chan))];  % UF-EN
%         end
%         for c = 1:size(y,1)
%             sy = smooth(y(c,:),length(y),'sgolay',9);
%             [namp,nlocs] = findpeaks(-sy(startpoint:end));
%             [pamp,plocs] = findpeaks(sy(startpoint:end));
%             peaktime(s,c,chan) = min([nlocs; plocs]) + startpoint;
%         end
%     end
end

figure
bar(mean(xcoef))

pvals = [];
for i = 1:size(combos,1)
    hold on; plot([i i],[mean(xcoef(:,i))+(std(xcoef(:,i))/sqrt(N)) mean(xcoef(:,i))-(std(xcoef(:,i))/sqrt(N))],'k');
    [h,p] = ttest(xcoef(:,i));
    pvals(i,1) = round(p,3);
end
set(gca,'xticklabels',{'N-F','E-U','EN-UN','EN-EF','EF-UF','UN-UF'})


figure
for chan = 1:size(peaktime,3)
    subplot(8,8,chan)
    y = squeeze(peaktime(:,:,chan));
    m = mean(y);
    sem = std(y)/sqrt(N);
    upper = m+sem;
    lower = m-sem;

    bar(m); hold on
    for c = 1:4
        plot([c c],[lower(c) upper(c)],'k');
    end
end

%% Group level

theseSubjects = setdiff(1:N,[3]); % subject 3 missed about half the trials
thisN = length(theseSubjects);
            
cfg = [];
cfg.method = 'template';
cfg.layout = 'biosemi64.lay';
neighbours = ft_prepare_neighbours(cfg);

if includeRT
    filename = 'individualStats_RT.mat';
else
    filename = 'individualStats.mat';
end

switch locking
    case 'stimulus'
        filename = ['SL_' filename];
    case 'response'
        filename = ['RL_' filename];
end

load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,filename)); % loads 'STATS' variable

% load('D:\bCFS_EEG_Reanalysis\data\Exp2\behav\stai.mat'); % loads 'stai' variable
% stai = table2array(struct2table(stai));
% stai = stai(theseSubjects,:);
% for i = 1:size(stai,2)
%     stai(:,i) = stai(:,i) - mean(stai(:,i));
% end
% stai = array2table(stai,'variablenames',{'full','state','trait'});

for st = 3%1:length(standardType)
    for b = 1%:length(baselineCorrect)
        stat = cell(1,3);
        for i = 1%:3 % 1: [UN-EF]-[UF-EN], 2: UN-EF, 3: UF-EN
            
            if strcmp(locking,'response')
                minrt = nan(N,3); % min time, max time, num samples
                for s = 1:N
                    minrt(s,:) = [STATS{s}{st,b,i}.time([1 end]) length(STATS{s}{st,b,i}.time)];
                end
                minrt = minrt(theseSubjects,:);
                minrt = max(minrt(:,1));
            end
            
            % Plot the prediction error for each emotion type
            if i==1
                
                C = 2;
                
                pdata = cell(N,2);
                for s = 1:N
                    for c = 1:C
                        thisstat = STATS{s}{st,b,c+1};
                        pdata{s,c}.dimord = thisstat.dimord;
                        pdata{s,c}.label = thisstat.label;
                        pdata{s,c}.time = thisstat.time;
                        pdata{s,c}.avg = zscore(thisstat.stat);
                    end
                end
                
                gpdata = cell(1,2);
                for c = 1:C
                    gpdata{c} = ft_timelockgrandaverage([],pdata{theseSubjects,c});
                end
                
                figure
                cfg = [];
                cfg.layout = 'biosemi64.lay';
                cfg.linecolor = cmap([2 4],:);
                cfg.linewidth = 1.4;
                ft_multiplotER(cfg,gpdata{:});
                
            end
            
            % Get the z-scored t maps for this statistic
            A = cell(N,1);
            for s = 1:N
                thisstat = STATS{s}{st,b,i};
                A{s,1}.dimord = thisstat.dimord;
                A{s,1}.label = thisstat.label;
                A{s,1}.time = thisstat.time;
                A{s,1}.avg = zscore(thisstat.stat);
                
                % For response-locking, make all subjects have the same length
                if strcmp(locking,'response')
                    thisx = findMin(minrt(1),A{s,1}.time):findMin(0,A{s,1}.time);
                    A{s,1}.time =  A{s,1}.time(thisx);
                    A{s,1}.avg = A{s,1}.avg(:,thisx);
                end
            end

            % Make dummy variable
            B = A;
            for s = 1:N
                B{s}.avg = zeros(size(B{s}.avg,1),size(B{s}.avg,2));
            end
            
            % Do t test
            cfg = [];
            cfg.neighbours = neighbours;
            cfg.channel ={'EEG','-M1','-M2'};

            cfg.method = 'montecarlo';
            cfg.tail = 0;
            cfg.alpha = .025;

            cfg.correctm = 'cluster';
            cfg.clusteralpha = .05;
            cfg.minnbchan = 0; % how many channels needed to form a cluster (minimum) - set to 0 to ignore this criterion
            cfg.numrandomization = 500;

            cfg.statistic = 'depsamplesT';
            cfg.design = [ones(1,thisN) ones(1,thisN)*2; 1:thisN 1:thisN];
            cfg.ivar = 1;
            cfg.uvar = 2;

            % ------------------------

            switch locking
                case 'stimulus'
                    cfg.latency = [0 3];
                case 'response'
                    cfg.latency = [minrt(1) 0];
            end

            cfg.avgoverchan = 'no';
            cfg.avgovertime = 'no';

            stat{i} = ft_timelockstatistics(cfg, A{theseSubjects}, B{theseSubjects});

            figure
            imagesc(1-stat{i}.prob);
            colormap('hot')
            xlabel('Time')
            ylabel('Channels')
            set(gca,'ytick',1:length(stat{i}.label));
            set(gca,'yticklabels',stat{i}.label);
            title(['Lowest p = ' num2str(min(stat{i}.prob(:)))])
            set(gcf,'position',[440 -141 484 957])
            drawnow;
            
            sigtime = {};
            sigchan = {};
            cc = 0;
            for posneg = 1:2
                
                if posneg==1 % positive clusters
                    thiscluster = stat{i}.posclusters;
                    thismap = stat{i}.posclusterslabelmat;
                elseif posneg==2 % negative clusters
                    thiscluster = stat{i}.negclusters;
                    thismap = stat{i}.negclusterslabelmat;
                end
                
                if ~isempty(thiscluster)
                    nclusters = find(extractfield(thiscluster,'prob') <= .05);
                    for n = 1:length(nclusters)
                        cc = cc + 1;
                        sigtime{cc} = stat{i}.time(sum(thismap == n) > 0);
                        sigchan{cc} = stat{i}.label(sum(thismap == n,2) > 0);
                    end
                end
            end
               
            figure
            for j = 1:length(sigtime)
                
                subplot(1,length(sigtime),j)
                
                cfg = [];
                cfg.parameter = 'stat';
                cfg.layout = 'biosemi64.lay';

                cfg.xlim = sigtime{j}([1 end]);

                cfg.highlight = 'on';
                cfg.highlightchannel = sigchan{j};
                cfg.highlightsymbol = '.';
                cfg.highlightcolor = [1 1 1];
                
                cfg.style = 'straight';
                cfg.markersymbol  = '.';
                cfg.markersize = 6;
                cfg.markercolor = [0 0 0];
                cfg.comment = 'no';

                ft_topoplotER(cfg,stat{i});
                colormap(colours(100,'viridis'))
                colorbar
                
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
%     tmp = load(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc0.mat']));
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

for s = 1:thisN
    
    subjslopes = slopes(slopes.subject==theseSubjects(s),:);
    
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
