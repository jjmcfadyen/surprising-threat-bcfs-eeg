clear all
clc

includeRT = false;
locking = 'stimulus'; % 'stimulus' or 'response'

%% Directories

dir_data = 'D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\preprocessed';
dir_imgs = 'D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\images';

subjects = [1:22 24:33]; % 1:33;
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
    for st = 3%1:length(standardType)
        
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
            switch locking
                case 'stimulus'
                    tmp = load(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc' num2str(baselineCorrect(b)) '.mat']));
                case 'response'
                    if ~baselineCorrect(b)
                        tmp = load(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc' num2str(baselineCorrect(b)) '.mat']));
                    else
                        % Use post-response 100 ms as baseline
                        tmp = load(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc0.mat']));
                        for trl = 1:length(tmp.thisD.trial)
                            thisrt = tmp.thisT.RT(trl);
                            thistime = tmp.thisD.time{trl};
                            thisbc = findMin(thisrt,thistime):findMin(thisrt+.1,thistime);
                            tmp.thisD.trial{trl} = tmp.thisD.trial{trl} - mean(tmp.thisD.trial{trl}(:,thisbc),2);
                        end
                    end
            end
            thisD = tmp.thisD;
            idx = tmp.idx;
            thisT = tmp.thisT;
            
            if size(thisT,1) ~= length(thisD.trial)
                error('Mismatch with behavioural file')
            end
            
            % Check for bad channels
            channels = find(~ismember(thisD.label,...
                {'M1','M2','Nz','SNz','LH','RH','LV','UV','Status'}));
            
            nTrls = length(thisD.trial);
            nChan = length(channels);
            
            d = nan(nChan,length(thisD.time{1}),nTrls);
            for trl = 1:nTrls
                d(:,:,trl) = thisD.trial{trl}(channels,:);
            end
            
            d = std(reshape(d,size(d,1),size(d,2)*size(d,3)),[],2); % get standard deviation of each channel
            
            outliers = gesd(d,.005); % relatively more strict alpha criterion for what makes a channel an outlier
            
            badchannels = thisD.label(channels(outliers));
            thisinterp{st,b} = badchannels;
            
            disp(['Interpolating: ' badchannels'])
            
            % Interpolate
            cfg = [];
            cfg.method = 'weighted';
            cfg.badchannel = badchannels;
            cfg.neighbours = neighbours;
            interp = ft_channelrepair(cfg,thisD);
            
            % Shift time axis for response-locked data
            if strcmp(locking,'response')
               
                nTrls = length(interp.trial);
                offset = nan(nTrls,1);
                for trl = 1:nTrls
                    offset(trl,1) = -findMin(thisT.RT(trl),interp.time{trl});
                end
                
                cfg = [];
                cfg.offset = offset;
                interp = ft_redefinetrial(cfg,interp);
                
                minmax = nan(nTrls,2);
                for trl = 1:nTrls
                    minmax(trl,:) = interp.time{trl}([1 end]); 
                end
                
                rmidx = minmax(:,2) < .1;  % remove responses too close to the end of the trial
                idx(rmidx,:) = 0;
                
                minmax = [max(minmax(idx,1)) min(minmax(idx,2))];
                
            end

            % Select data
            data = cell(1,4);
            for c = 1:4
                cfg = [];
                cfg.trials = find(idx & thisT.Condition==c);
                data{c} = ft_selectdata(cfg,interp);
            end
            
            % If response-locked, make all trials the same length
            if strcmp(locking,'response')
                cfg = [];
                cfg.toilim = minmax;
                for c = 1:4
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
                    cfg.latency = [minmax(1) 0];
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
                            A.trial{trl} = A.trial{trl} - avEF.avg; % minus average EF
                        end
                        
                        B = data{4}; % UF
                        for trl = 1:length(B.trial)
                            B.trial{trl} = B.trial{trl} - avEN.avg; % minus average EN
                        end

                    elseif i==2 % UN vs EF

                        A = data{2}; % UN
                        B = data{3}; % EF

                    elseif i==3 % UF vs EN

                        A = data{4}; % UF
                        B = data{1}; % EN

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

%% Group level

theseSubjects = setdiff(1:N,[3 5]); % subject 3 missed about half the trials, 11 & 24 are outlier amplitude (noisy)
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

load(fullfile('D:\bCFS_EEG_Reanalysis\results\sensorspace',filterTag,filename));

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
                minrt = max(minrt);
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
            
            % Plot the z-scored t maps for this statistic
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
            
%             if ~includeRT
                cfg.channel = {'EEG','-M1','-M2'};
%             else
%                 cfg.channel = {'F1','F3','F5','FC3','C1','C3','C5','T7','TP7','CP5','CP3','P3','P5','P7','P9','PO7','O1','Iz','Oz','Fp2','AF8','AF4','F4','F8','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP4','CP2','P2','P8','P10','PO8','O2'};
%             end

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
                cfg.highlightcolor = [1 1 1];

                ft_topoplotER(cfg,stat{i});
                colormap(colours(100,'viridis'))
                colorbar
                
            end
        end
    end
end
