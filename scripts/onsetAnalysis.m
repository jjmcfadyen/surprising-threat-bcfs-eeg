% Analyse the onset times
clear all
clc

%% Directories

addpath('D:\Toolboxes\fieldtrip-master')
ft_defaults;

dir_data = 'D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\preprocessed\unfiltered';

%% Parameters

driftwindow = [0.5 0.1];
nIterations = 0; % if 0, then don't do randomisation of response times

subjects = [1:2 4:22 24:33]; % 1:33;
N = length(subjects);

cmap = [77, 242, 255;
        98, 174, 255;
        255, 180, 29;
        255, 68, 142]/255;

cfg = [];
cfg.method = 'template';
cfg.layout = 'biosemi64.lay';
neighbours = ft_prepare_neighbours(cfg);

standardType        = {'all','first', 'last'};
baselineCorrect     = [true, false];
st = 3; % last standards
b = 1;  % baseline corrected epochs
datatype = 'orig';

%% Estimate onsets for each trial

Y = [];
rts = nan(N,4);
missingdata = nan(N,64,2,nIterations+1);
onsets = nan(N,4,64,2);
rlonsets = onsets;
nullonsets = onsets;
nullrlonsets = onsets;
alltrialonsets = cell(N,4,nIterations+1);
allrltrialonsets = alltrialonsets;
rtcorr = nan(N,64,2,nIterations+1);
rlrtcorr = nan(N,64,2,nIterations+1);
onsetdata = [];
rlonsetdata = [];
trialdrifts = nan(N,4,64);
alltrialdrifts = cell(N,4);
for s = 1:N

    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end
    
    disp(['Reading in data for ' schar '...'])
    
    % Load data
    [SL,RL,T,idx] = select_data(fullfile(dir_data,'newepoched',[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),neighbours,standardType{st},datatype);
    
    fidx = find(idx);
    nTrls = length(fidx);
    nChan = length(SL.label);

    for c = 1:4
        rts(s,c) = nanmean(T.RT(idx & T.Condition==c));
    end

    % Get data matrix & rate of change
    X = SL.time{1};
    xidx = X<=3;

    d = [];
    for trl = 1:nTrls
        if trl==1
            roc = nan(nTrls,nChan,sum(xidx));
        end
        for chan = 1:nChan
            tmp = SL.trial{fidx(trl)}(chan,xidx);
            tmp = movmean(tmp,10); % smooth
            d(trl,chan,:) = tmp;
        end
    end

    % See when each trial & each channel exceeds a certain threshold
    tic
    [trlonsets,rltrlonsets] = estimateOnsets(d,X(xidx),T.RT(idx),T.Condition(idx),nIterations);
    toc

    orig = [];
    orig.trlonsets = squeeze(trlonsets(:,:,:,1));
    orig.rltrlonsets = squeeze(rltrlonsets(:,:,:,1));

    % correlate each iteration with RT
    for i = 1:nIterations+1
        for chan = 1:nChan
            for j = 1:2

                thisx = T.RT(idx);
                thisy = squeeze(trlonsets(:,chan,j,i));
                nanidx = isnan(thisy);
                thisx = thisx(~nanidx);
                thisy = thisy(~nanidx);
                if isempty(thisy)
                    rtcorr(s,chan,j,i) = NaN;
                else
                    rtcorr(s,chan,j,i) = corr(thisx,thisy);
                end

                thisx = T.RT(idx);
                thisy = squeeze(rltrlonsets(:,chan,j,i));
                nanidx = isnan(thisy);
                thisx = thisx(~nanidx);
                thisy = thisy(~nanidx);
                if isempty(thisy)
                    rlrtcorr(s,chan,j,i) = NaN;
                else
                    rlrtcorr(s,chan,j,i) = corr(thisx,thisy);
                end
            end
        end
    end

    % log
    missingdata(s,:,:,i) = mean(isnan(trlonsets(:,:,:,1)));
    for c = 1:4
        onsets(s,c,:,:) = nanmean(trlonsets(T.Condition(idx)==c,:,:,1));
        rlonsets(s,c,:,:) = nanmean(rltrlonsets(T.Condition(idx)==c,:,:,1));
        if nIterations>0
            nullonsets(s,c,:,:) = nanmean(nanmean(trlonsets(T.Condition(idx)==c,:,:,2:end),4));
            nullrlonsets(s,c,:,:) = nanmean(nanmean(rltrlonsets(T.Condition(idx)==c,:,:,2:end),4));
        end
        for i = 1:nIterations+1
            alltrialonsets{s,c,i} = squeeze(trlonsets(T.Condition(idx)==c,:,:,i));
            allrltrialonsets{s,c,i} = squeeze(rltrlonsets(T.Condition(idx)==c,:,:,i));
        end
    end

    % Gather data to plot, organised by onset time
    thisx = X(xidx);
    for j = 1:2 % first onset, last onset
        for chan = 1:nChan
            q = quantile(orig.trlonsets(:,chan,j),4);
            for i = 1:3
    
                % get index of trials with this onset time (fast, medium, or slow)
                qidx = orig.trlonsets(:,chan,j)>q(i) & orig.trlonsets(:,chan,j)<q(i+1);
                fqidx = find(qidx);
    
                % select stimulus-locked data for this onset group and average it
                onsetdata(s,i,chan,:,j) = squeeze(mean(d(qidx,chan,:)));
    
                % log average onset times for this channel/quantile
                rtonsetdata(s,i,chan,1,j) = mean(orig.trlonsets(qidx,chan,j));
    
            end
        end
    end

    % Gather response-locked data to plot, organised by onset time
    thisx = X(xidx);
    theserts = T.RT(idx);
    for otype = 1:2 % first onset, last onset
        for chan = 1:nChan
            q = fliplr(quantile(orig.rltrlonsets(:,chan,otype),4)); % fastest to slowest
            for i = 1:3
    
                % get index of trials with this onset time (fast, medium, or slow)
                qidx = orig.rltrlonsets(:,chan,otype)<q(i) & orig.rltrlonsets(:,chan,otype)>q(i+1);
                fqidx = find(qidx);
    
                % select response-locked data
                tmp = nan(sum(qidx),size(d,3));
                for j = 1:length(fqidx)
                    thisrt = theserts(fqidx(j));
                    thisrltrial = squeeze(d(fqidx(j),chan,thisx<=(thisrt+0.1))); % get the data from first time point to response (plus 100 ms after)
                    tmp(j,size(tmp,2)-length(thisrltrial)+1:end) = thisrltrial; % insert into 'tmp' matrix, with leading NaNs
                end
                thisrlmean = nan(1,size(tmp,2));
                thisrlmean(:,mean(isnan(tmp))<0.05) = nanmean(tmp(:,mean(isnan(tmp))<0.05)); % save mean for time points where there is at least 95% of the data
                rlonsetdata(s,i,chan,:,otype) = thisrlmean;
    
                % log average onset times for this channel/quantile
                rtonsetdata(s,i,chan,2,otype) = mean(orig.rltrlonsets(qidx,chan,otype));
    
            end
        end
    end

    % Get drift rate per trial using predefined window
    v = nan(nTrls,nChan);
    for trl = 1:nTrls
        thisrt = T.RT(fidx(trl));
        thisx = SL.time{fidx(trl)}(xidx);
        thiswindow = findMin(thisrt-driftwindow(1),thisx):findMin(thisrt-driftwindow(2),thisx); % -800 to -200 before response onset
        for chan = 1:nChan
            thisy = squeeze(d(trl,chan,thiswindow));
            coeff = polyfit(thisx(thiswindow),thisy,1);
            v(trl,chan) = abs(coeff(1));
        end
    end

    for c = 1:4
        trialdrifts(s,c,:) = mean(v(T.Condition(idx)==c,:));
        alltrialdrifts{s,c} = v(T.Condition(idx)==c,:);
    end

    % Save average rate of change waves
    for chan = 1:nChan
        for c = 1:4
            Y(s,c,chan,:) = squeeze(mean(d(T.Condition(fidx)==c,chan,:)));
            ROC(s,c,chan,:) = squeeze(nanmean(roc(T.Condition(fidx)==c,chan,:)));
        end
    end
end

labels = SL.label;
x = X(xidx);

save('erponsets.mat',...
    'onsets','rlonsets','alltrialonsets','allrltrialonsets','rts','missingdata','onsetdata','rlonsetdata','rtonsetdata','driftwindow','labels','x');

% plot histograms
% get all trial onset times and convert to percentage
y = [];
ypcnt = [];
for s = 1:N

    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end
    [~,~,T,idx] = select_data(fullfile(dir_data,'newepoched',[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),neighbours,standardType{st},datatype);
    
    tmp = alltrialonsets{s,1};
    tmprt = T.RT(T.Condition==1 & idx);
    for c = 2:4
        cc = [size(tmp,1)+1 size(tmp,1)+size(alltrialonsets{s,c},1)];
        tmp(cc(1):cc(2),:,:) = alltrialonsets{s,c};
        tmprt = [tmprt; T.RT(T.Condition==c & idx)];
    end
    tmp = squeeze(mean(tmp(:,ismember(labels,{'CPz','Cz','Pz','PO4','PO8','O2','PO3','PO7','O1','Oz','POz'}),:),2));

    y = [y; tmp tmprt];
    ypcnt = [ypcnt; tmp ./ tmprt];
end

cmap = [243, 140, 255 ;
        66, 0, 213]/255;

figure
subplot(2,1,1)
for c = 1:2

    [bincounts,binedges] = histcounts(y(:,c),'binwidth',0.1);
    bincenters = binedges(1:end-1)+diff(binedges)/2;
    
    bar(bincenters,bincounts,1,'facecolor',cmap(c,:),'facealpha',0.5); hold on

    pd = fitdist(y(:,c),'normal');
    q = icdf(pd,[0.0013499 0.99865]); % three-sigma range for normal distribution
    thisx = linspace(q(1),q(2));

    binwidth = binedges(2)-binedges(1); % Finds the width of each bin
    area = size(y,1) * binwidth;
    thisfit = area * pdf(pd,thisx);

    plot(thisx,thisfit,'color',cmap(c,:),'linewidth',1.3); hold on
end
set(gca,'ticklength',[0 0])
xlim([0 3])

subplot(2,1,2)
for c = 1:2

    [bincounts,binedges] = histcounts(ypcnt(:,c),'binwidth',0.05);
    bincenters = binedges(1:end-1)+diff(binedges)/2;
    
    bar(bincenters,bincounts,1,'facecolor',cmap(c,:),'facealpha',0.5); hold on

    pd = fitdist(ypcnt(:,c),'normal');
    q = icdf(pd,[0.0013499 0.99865]); % three-sigma range for normal distribution
    thisx = linspace(q(1),q(2));

    binwidth = binedges(2)-binedges(1); % Finds the width of each bin
    area = size(ypcnt,1) * binwidth;
    thisfit = area * pdf(pd,thisx);

    plot(thisx,thisfit,'color',cmap(c,:),'linewidth',1.3); hold on
end
set(gca,'ticklength',[0 0])
xlim([0 1])

%% Plot

chanselection = [];
chanselection{1} = {'PO4','PO8','O2','PO3','PO7','O1','Oz','POz'};
chanselection{2} = {'CPz','Cz','Pz'};

cmap = [243, 140, 255 ;
        178, 71, 255 ;
        66, 0, 213]/255;

for j = 1:2
    figure
    for chan = 1:length(chanselection) % specific channel groupings
    
        subplot(1,length(chanselection),chan)
    
        for i = 1:3
    
            y = squeeze(mean(onsetdata(:,i,ismember(labels,chanselection{chan}),:,j),3));
            m = mean(y);
            sem = std(y)/sqrt(size(y,1));
            upper = m+sem;
            lower = m-sem;
    
%             patch([x fliplr(x)],[upper fliplr(lower)],cmap(i,:),'facealpha',.2,'edgealpha',0); hold on
            plot(x,m,'color',cmap(i,:),'linewidth',1.3); hold on
        end
        
        ax = gca;
        for i = 1:3
            plot(repmat(mean(rtonsetdata(:,i,chan,1,j)),2,1),ax.YLim,'color',cmap(i,:),'linewidth',1.2,'linestyle','--'); hold on
        end
        plot(x([1 end]),[0 0],'k:'); hold on
        title(strjoin(chanselection{chan}))
        xlim(x([1 end]))
    
    end
    if j==1
        sgtitle('SL early onset times (fast, medium, slow)')
    elseif j==2
        sgtitle('SL late onset times (fast, medium, slow)')
    end
end

% -- response-locked
thisx = x - 2.9;
for j = 1:2
    figure
    for chan = 1:length(chanselection) % specific channel groupings
    
        subplot(1,length(chanselection),chan)
    
        for i = 1:3
    
            y = squeeze(nanmean(rlonsetdata(:,i,ismember(labels,chanselection{chan}),:,j),3));
            m = nanmean(y);
            sem = nanstd(y)/sqrt(size(y,1));
            upper = m+sem;
            lower = m-sem;
    
            nanidx = isnan(upper) | isnan(lower);
%             patch([thisx(~nanidx) fliplr(thisx(~nanidx))],[upper(~nanidx) fliplr(lower(~nanidx))],cmap(i,:),'facealpha',.2,'edgealpha',0); hold on
            plot(thisx,m,'color',cmap(i,:),'linewidth',1.3); hold on
        end
        
        ax = gca;
        for i = 1:3
            plot(repmat(nanmean(rtonsetdata(:,i,chan,2,j)),2,1),ax.YLim,'color',cmap(i,:),'linewidth',1.2,'linestyle','--'); hold on
        end
        tmpx = ~isnan(squeeze(nanmean(nanmean(rlonsetdata(:,:,chan,:,j)),2)));
        tmpx = thisx([find(tmpx,1,'first') find(tmpx,1,'last')]);
        plot(tmpx,[0 0],'k:'); hold on
        plot([0 0],ax.YLim,'k:'); hold on
        title(strjoin(chanselection{chan}))
%         xlim([-1.5 0.1])
    
    end
    
    if j==1
        sgtitle('RL early onset times (fast, medium, slow)')
    elseif j==2
        sgtitle('RL late onset times (fast, medium, slow)')
    end
end

%% Save trial onsets to long table for analysis in R

chanselection = [];
chanselection{1} = {'PO4','PO8','O2','PO3','PO7','O1','Oz','POz'};
chanselection{2} = {'CPz','Cz','Pz'};

changroupnames = {'PO','CP'};

ddm = readtable('D:\bCFS_EEG_Reanalysis\results\ddmparams_nozscore.csv'); 
ddm.Subject(ddm.Experiment==2) = ddm.Subject(ddm.Experiment==2)-31; % 31 subjects in experiment 1
ddm = ddm(ddm.Experiment==2,:);
ddm = ddm(ismember(ddm.Subject,subjects),:);

load('D:\bCFS_EEG_Reanalysis\data\Exp2\behav\stai.mat'); % loads 'stai' variable
stai = table2array(struct2table(stai));
stai = stai(subjects,:);
stai = array2table(stai,'variablenames',{'full','state','trait'});

L = [];
for s = 1:N

    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end

    [~,~,T,idx] = select_data(fullfile(dir_data,'newepoched',[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),neighbours,standardType{st},datatype);

    for c = 1:4 % trials batched into conditions

        nTrls = size(alltrialonsets{s,c},1);

        for chan = 1:length(chanselection) % electrode type
            for o = 1:3 % change point type (early, late, or difference)
    
                thistable = array2table([repmat(s,nTrls,1) repmat(c,nTrls,1)],'variablenames',{'Subject','Condition'});
                thistable.Trial = find(T.Condition==c & idx);

                if c <= 2
                    thistable.Emotion = repmat({'Neutral'},nTrls,1);
                else
                    thistable.Emotion = repmat({'Fearful'},nTrls,1);
                end
                if c ==1 || c==3
                    thistable.Expectation = repmat({'Expected'},nTrls,1);
                else
                    thistable.Expectation = repmat({'Unexpected'},nTrls,1);
                end

                if chan==1
                    thistable.Electrode = repmat({'parietaloccipital'},nTrls,1);
                elseif chan==2
                    thistable.Electrode = repmat({'centralparietal'},nTrls,1);
                end

                if o==1
                    thistable.onsetType = repmat({'early'},nTrls,1);
                elseif o==2
                    thistable.onsetType = repmat({'late'},nTrls,1);
                elseif o==3
                    thistable.onsetType = repmat({'diff'},nTrls,1);
                end
        
                if o<=2
                    thistable.onsetSL = mean(alltrialonsets{s,c}(:,ismember(labels,chanselection{chan}),o),2);
                    thistable.onsetRL = mean(allrltrialonsets{s,c}(:,ismember(labels,chanselection{chan}),o),2);
                elseif o==3
                    thistable.onsetSL = mean(alltrialonsets{s,c}(:,ismember(labels,chanselection{chan}),2),2) - mean(alltrialonsets{s,c}(:,ismember(labels,chanselection{chan}),1),2);
                    thistable.onsetRL = mean(allrltrialonsets{s,c}(:,ismember(labels,chanselection{chan}),2),2) - mean(allrltrialonsets{s,c}(:,ismember(labels,chanselection{chan}),1),2);
                end
                thistable.slope = mean(alltrialdrifts{s,c}(:,ismember(labels,chanselection{chan})),2);
        
                thistable.RT = T.RT(idx & T.Condition==c);
        
                thistable.nondecision = nan(nTrls,1);
                thistable.drift = nan(nTrls,1);
                thistable.boundary = nan(nTrls,1);
                if any(ismember(ddm.Subject,subjects(s)))
                    thistable.nondecision = repmat(table2array(ddm(ismember(ddm.Subject,subjects(s)),ismember(ddm.Properties.VariableNames,['nondecision_C' num2str(c)]))),nTrls,1);
                    thistable.drift = repmat(table2array(ddm(ismember(ddm.Subject,subjects(s)),ismember(ddm.Properties.VariableNames,['drift_C' num2str(c)]))),nTrls,1);
                    thistable.boundary = repmat(table2array(ddm(ismember(ddm.Subject,subjects(s)),ismember(ddm.Properties.VariableNames,['boundary_C' num2str(c)]))),nTrls,1);
                end
        
                thistable.Anxiety = repmat(stai.trait(s),nTrls,1);
        
                L = [L; thistable];
            end
        end
    end
end

writetable(L,'D:\bCFS_EEG_Reanalysis\results\erponsets_alltrials.csv');

%%
% Onsets & drifts
ddm = readtable('D:\bCFS_EEG_Reanalysis\results\ddmparams_nozscore.csv'); 
ddm.Subject(ddm.Experiment==2) = ddm.Subject(ddm.Experiment==2)-31; % 31 subjects in experiment 1
ddm = ddm(ddm.Experiment==2,:);
ddm = ddm(ismember(ddm.Subject,subjects),:);

cmap = [0 232 255;
        51 110 255;
        255 186 48;
        255 0 0]/255;

for j = 1%:2 % stimulus-locked, response-locked
    if j==1
        locktype = 'SL';
    elseif j==2
        locktype = 'RL';
    end
    for i = 1:3 % channel group

        if i==1
            chanselection = {'PO4','PO8','O2','PO3','PO7','O1','Oz','POz'};
        elseif i==2
            chanselection = {'CPz','Cz','Pz'};
        elseif i==3
            chanselection = {'PO4','PO8','O2','PO3','PO7','O1','Oz','POz','CPz','Cz','Pz'};
        end

        tmp = squeeze(mean(rtcorr(:,ismember(labels,chanselection),:),2));
        [~,p,~,tstat] = ttest(tmp(:,1),tmp(:,2));
        disp(['SL Correlations of each onset type with RT: EARLY: ' num2str(round(mean(tmp(:,1)),3)) ', LATE: ' num2str(round(mean(tmp(:,2)),3)),...
            ', t = ' num2str(round(tstat.tstat,3)) ', p = ' num2str(round(p,3))])

        tmp = squeeze(mean(rlrtcorr(:,ismember(labels,chanselection),:),2));
        [~,p,~,tstat] = ttest(tmp(:,1),tmp(:,2));
        disp(['RL Correlations of each onset type with RT: EARLY: ' num2str(round(mean(tmp(:,1)),3)) ', LATE: ' num2str(round(mean(tmp(:,2)),3)),...
            ', t = ' num2str(round(tstat.tstat,3)) ', p = ' num2str(round(p,3))])

        for otype = 1:2 % early, late
        
            if otype==1
                otitle = 'early';
            elseif otype==2
                otitle = 'late';
            end

            thisy = [];
            if j==1 % stimulus-locked onsets
                thisy = squeeze(mean(onsets(:,:,ismember(labels,chanselection),otype),3));
            elseif j==2 % response-locked onsets
                thisy = squeeze(mean(rlonsets(:,:,ismember(labels,chanselection),otype),3));
            end
            thisy = zscore(thisy')';
        
            figure
            m = mean(thisy);
            sem = std(thisy)/sqrt(size(thisy,1));
            upper = m+sem;
            lower = m-sem;
            for c = 1:4
        
                q = quantile(thisy(:,c),[.25 .75]);
                patch([c-.15 c-.15 c+.15 c+.15],[q(1) q(2) q(2) q(1)],cmap(c,:),'facealpha',.25,'edgecolor',cmap(c,:)); hold on
                plot([c c],[q(2) max(thisy(:,c))],'color',cmap(c,:),'linewidth',1.2); hold on
                plot([c c],[q(1) min(thisy(:,c))],'color',cmap(c,:),'linewidth',1.2); hold on
        
                plot([c c],[lower(c) upper(c)],'k','linewidth',1.3); hold on
                scatter(c,m(c),70,'markerfacecolor',cmap(c,:),'markeredgecolor','k'); hold on
            end
            [~,p] = ttest(thisy(:,2)-thisy(:,1),thisy(:,4)-thisy(:,3));
            title([locktype ', ' otitle ', ' strjoin(chanselection,' ') ', p = ' num2str(round(p,3))])
            xlim([0 5])
            set(gca,'ticklength',[0 0])
%             set(gcf,'position',[715 332 296 336])
        
        %     % save
        %     writetable(array2table(thisy,'variablenames',{'EN','UN','EF','UF'}),fullfile('D:\bCFS_EEG_Reanalysis\results',['erponsets_changroup' num2str(i) '.csv']))
        
            ST = array2table([[1:N]' thisy],'variablenames',{'subject','neutral_expected','neutral_unexpected','fearful_expected','fearful_unexpected'});
            ST = stack(ST,{'neutral_expected','neutral_unexpected','fearful_expected','fearful_unexpected'});
            ST.Properties.VariableNames = {'subject','condition','y'};
            ST.emotion = double(ST.condition=='fearful_expected' | ST.condition=='fearful_unexpected');
            ST.expectation = double(ST.condition=='neutral_unexpected' | ST.condition=='fearful_unexpected');
            stats = rm_anova2(ST.y,ST.subject,ST.emotion,ST.expectation,{'Emotion','Expectation'});
            stats = array2table(stats(2:end,:),'variablenames',stats(1,:))
        
            [~,p] = ttest(thisy(:,1),thisy(:,2));
            disp(['EN vs UN: p = ' num2str(round(p,3))])
        
            [~,p] = ttest(thisy(:,3),thisy(:,4));
            disp(['EF vs UF: p = ' num2str(round(p,3))])

            % correlate with DDM parameters
            for k = 1:3
        
                figure
        
                if k==1
                    thisparam = 'nondecision';
                elseif k==2
                    thisparam = 'drift';
                elseif k==3
                    thisparam = 'boundary';
                elseif k==4
                    thisparam = 'rt';
                end

                if k<=3
                    if j==1 % stimulus-locked onsets
                        thiseeg = squeeze(mean(onsets(ismember(subjects,ddm.Subject),:,ismember(labels,chanselection),otype),3));
                    elseif j==2 % response-locked onsets
                        thiseeg = squeeze(mean(rlonsets(ismember(subjects,ddm.Subject),:,ismember(labels,chanselection),otype),3));
                    end
                    thisddm = table2array(ddm(:,contains(ddm.Properties.VariableNames,thisparam)));
                elseif k==4
                    if j==1 % stimulus-locked onsets
                        thiseeg = squeeze(mean(onsets(:,:,ismember(labels,chanselection),otype),3));
                    elseif j==2 % response-locked onsets
                        thiseeg = squeeze(mean(rlonsets(:,:,ismember(labels,chanselection),otype),3));
                    end
                    thisddm = rts;
                end
                thiseeg = zscore(thiseeg')';
                thisddm = zscore(thisddm')';

                scatter(thiseeg(:),thisddm(:),'markerfacecolor','k','markeredgecolor','none','markerfacealpha',.5); hold on
                coeff = polyfit(thiseeg(:),thisddm(:),1);
                thisline = polyval(coeff,thiseeg(:));
                plot(thiseeg(:),thisline,'k'); hold on
                [r,p] = corr(thiseeg(:),thisddm(:));
                title([locktype ', ' otitle ', ' strjoin(chanselection,' ') ', ' thisparam ': r=' num2str(round(r,3)) ', p=' num2str(round(p,3))])
            end
        end

%         % correlation with RT
%         thisy = squeeze(mean(rtcorr(:,ismember(labels,chanselection),:),2));
%         m = mean(thisy);
%         sem = std(thisy)/sqrt(size(thisy,1));
%         upper = m+sem;
%         lower = m-sem;
%         figure
%         bar(m); hold on
%         for c = 1:2
%             plot([c c],[upper(c) lower(c)],'k'); hold on
%         end
%         set(gca,'xticklabels',{'early','late'})
%         [~,p,~,tstats] = ttest(thisy(:,1),thisy(:,2));
%         title([locktype ', ' strjoin(chanselection,' ') ', onset-RT correlation: p = ' num2str(p)])

    end
end

%% Cluster-based permutation 

x = X(xidx);
thisd = ROC; % Y or ROC

stat = [];
for i = 1:5 % emotion, expectation, interaction, neut pred, fear pred

    A = cell(1,N);
    B = cell(1,N);
    for s = 1:N

        A{s} = ft_timelockanalysis([],SL);
        A{s}.var = nan(size(thisd,3),size(thisd,4));
        A{s}.dof = nan(size(thisd,3),size(thisd,4));
        A{s}.time = x;
        B{s} = A{s};

        if i==1
            A{s}.avg = squeeze(mean(thisd(s,1:2,:,:),2));
            B{s}.avg = squeeze(mean(thisd(s,3:4,:,:),2));
        elseif i==2
            A{s}.avg = squeeze(mean(thisd(s,[1 3],:,:),2));
            B{s}.avg = squeeze(mean(thisd(s,[2 4],:,:),2));
        elseif i==3
            A{s}.avg = squeeze(thisd(s,2,:,:))-squeeze(thisd(s,1,:,:));
            B{s}.avg = squeeze(thisd(s,4,:,:))-squeeze(thisd(s,3,:,:));
        elseif i==4
            A{s}.avg = squeeze(thisd(s,2,:,:));
            B{s}.avg = squeeze(thisd(s,1,:,:));
        elseif i==5
            A{s}.avg = squeeze(thisd(s,4,:,:));
            B{s}.avg = squeeze(thisd(s,3,:,:));
        end
    end

    % Do t test
    cfg = [];
    cfg.neighbours = neighbours;
    cfg.channel ={'EEG','-M1','-M2'};
%     cfg.channel = {'C1','C2','Cz','CPz','CP1','CP2','Pz','P1','P2'};

    cfg.method = 'montecarlo';
    cfg.tail = 0;
    cfg.alpha = .025;

    cfg.correctm = 'cluster';
    cfg.clusteralpha = .05;
    cfg.minnbchan = 2;
    cfg.numrandomization = 500;

    cfg.statistic = 'depsamplesT'; % 'indepsamplesregrT'
    cfg.design = [ones(1,N) ones(1,N)*2; 1:N 1:N];
    cfg.ivar = 1;
    cfg.uvar = 2;

    cfg.latency = [0 3];

    cfg.avgoverchan = 'no';
    cfg.avgovertime = 'no';

    stat = ft_timelockstatistics(cfg, A{:}, B{:});

    figure
    imagesc(1-stat.prob);
    colormap('hot')
    xlabel('Time')
    ylabel('Channels')
    set(gca,'ytick',1:length(stat.label));
    set(gca,'yticklabels',stat.label);
    title(['Lowest p = ' num2str(min(stat.prob(:)))])
    set(gcf,'position',[440 -141 484 957])
    drawnow;
end

%% Cross-correlate

thiscorr = [];
pvals = [];
for chan = 1:nChan
    for s = 1:N
        thisy = squeeze(Y(s,:,chan,:));
        if all(isnan(thisy(:)))
            thiscorr(s,chan) = NaN;
        else
            
            [C,lags] = xcorr(thisy(1,:),thisy(2,:));
            thiscorr(s,chan,1) = lags(C==max(C));

            [C,lags] = xcorr(thisy(3,:),thisy(4,:));
            thiscorr(s,chan,2) = lags(C==max(C));

        end
    end
    [~,p] = ttest(squeeze(thiscorr(:,chan,1)),squeeze(thiscorr(:,chan,2)));
    pvals(chan) = p;
end

[sorted,sortidx] = sort(pvals);
bestChan = A{1}.label(sortidx(sorted < .05));

%% Cluster-based permutation with Ter as DDM regressor

ddm = readtable('D:\bCFS_EEG_Reanalysis\results\ddmparams_nozscore.csv');
ddm.Subject(ddm.Experiment==2) = ddm.Subject(ddm.Experiment==2)-31; % 31 subjects in experiment 1
ddm = ddm(ddm.Experiment==2,:);
ddm = ddm(ismember(ddm.Subject,subjects),:);

ter = [ddm.nondecision_C1 ddm.nondecision_C2 ddm.nondecision_C3 ddm.nondecision_C4];
v = [ddm.drift_C1 ddm.drift_C2 ddm.drift_C3 ddm.drift_C4];

for i = 1:7 % emotion, expectation, interaction, neut pred, fear pred, pool conditions, average across conditions

    thisddm = ter;
    thisd = Y; % Y or ROC

    A = cell(1,N);
    for s = 1:N

        A{s} = ft_timelockanalysis([],SL);
        A{s}.var = nan(size(thisd,3),size(thisd,4));
        A{s}.dof = nan(size(thisd,3),size(thisd,4));
        A{s}.time = x;

        if i==1
            A{s}.avg = squeeze(mean(ROthisdC(s,1:2,:,:),2)) - squeeze(mean(thisd(s,3:4,:,:),2));
        elseif i==2
            A{s}.avg = squeeze(mean(thisd(s,[1 3],:,:),2)) - squeeze(mean(thisd(s,[2 4],:,:),2));
        elseif i==3
            A{s}.avg = (squeeze(thisd(s,2,:,:)) - squeeze(thisd(s,1,:,:))) - (squeeze(thisd(s,4,:,:)) - squeeze(thisd(s,3,:,:)));
        elseif i==4
            A{s}.avg = (squeeze(thisd(s,2,:,:)) - squeeze(thisd(s,1,:,:)));
        elseif i==5
            A{s}.avg = (squeeze(thisd(s,4,:,:)) - squeeze(thisd(s,3,:,:)));
        elseif i==6
            tmp = A{s};
            for c = 1:4
                A{c,s} = tmp;
                A{c,s}.avg = squeeze(thisd(s,c,:,:));
            end
        elseif i==7
            A{s}.avg = squeeze(mean(thisd(s,:,:,:),2));
        end
    end
    A = A(:,ismember(subjects,ddm.Subject))';
    A = A(:);

    B = nan(1,size(thisddm,1));
    for s = 1:size(thisddm,1)
        if i==1
            B(s) = mean(thisddm(s,[1 2]),2) - mean(thisddm(s,[3 4]),2);
        elseif i==2
            B(s) = mean(thisddm(s,[1 3]),2) - mean(thisddm(s,[2 4]),2);
        elseif i==3
            B(s) = (mean(thisddm(s,2),2)-mean(thisddm(s,1),2)) - (mean(thisddm(s,4),2)-mean(thisddm(s,3),2));
        elseif i==4
            B(s) = mean(thisddm(s,2),2)-mean(thisddm(s,1),2);
        elseif i==5
            B(s) = mean(thisddm(s,4),2)-mean(thisddm(s,3),2);
        elseif i==6
            for c = 1:4
                B(c,s) = thisddm(s,c);
            end
        elseif i==7
            B(s) = mean(thisddm(s,:));
        end
    end
    B = B(:);

    % Do t test
    cfg = [];
    cfg.neighbours = neighbours;
    cfg.channel ={'EEG','-M1','-M2'};
%     cfg.channel = {'C1','C2','Cz','CPz','CP1','CP2','Pz','P1','P2'};

    cfg.method = 'montecarlo';
    cfg.tail = 0;
    cfg.alpha = .025;

    cfg.correctm = 'cluster';
    cfg.clusteralpha = .05;
    cfg.minnbchan = 2;
    cfg.numrandomization = 500;

    cfg.statistic = 'ft_statfun_correlationT'; % 'indepsamplesregrT'
    cfg.design = B;
    cfg.ivar = 1;

    cfg.latency = [0 3];

    cfg.avgoverchan = 'no';
    cfg.avgovertime = 'no';

    stat = ft_timelockstatistics(cfg, A{:});

    figure
    imagesc(1-stat.prob);
    colormap('hot')
    xlabel('Time')
    ylabel('Channels')
    set(gca,'ytick',1:length(stat.label));
    set(gca,'yticklabels',stat.label);
    title(['Lowest p = ' num2str(min(stat.prob(:)))])
    set(gcf,'position',[440 -141 484 957])
    drawnow;

end

%% Cluster-based permutation with anxiety as regressor

load('D:\bCFS_EEG_Reanalysis\data\Exp2\behav\stai.mat'); % loads 'stai' variable
stai = table2array(struct2table(stai));
stai = stai(subjects,:);
for i = 1:size(stai,2)
    stai(:,i) = stai(:,i) - mean(stai(:,i));
end
stai = array2table(stai,'variablenames',{'full','state','trait'});

thisd = ROC;

for i = 1:5 % emotion, expectation, interaction, neut pred, fear pred

    A = cell(1,N);
    for s = 1:N

        A{s} = ft_timelockanalysis([],SL);
        A{s}.var = nan(size(thisd,3),size(thisd,4));
        A{s}.dof = nan(size(thisd,3),size(thisd,4));
        A{s}.time = x;

        if i==1
            A{s}.avg = squeeze(mean(thisd(s,1:2,:,:),2)) - squeeze(mean(thisd(s,3:4,:,:),2));
        elseif i==2
            A{s}.avg = squeeze(mean(thisd(s,[1 3],:,:),2)) - squeeze(mean(thisd(s,[2 4],:,:),2));
        elseif i==3
            A{s}.avg = (squeeze(thisd(s,2,:,:)) - squeeze(thisd(s,1,:,:))) - (squeeze(thisd(s,4,:,:)) - squeeze(thisd(s,3,:,:)));
        elseif i==4
            A{s}.avg = (squeeze(thisd(s,2,:,:)) - squeeze(thisd(s,1,:,:)));
        elseif i==5
            A{s}.avg = (squeeze(thisd(s,4,:,:)) - squeeze(thisd(s,3,:,:)));
        elseif i==6
            tmp = A{s};
            for c = 1:4
                A{c,s} = tmp;
                A{c,s}.avg = squeeze(thisd(s,c,:,:));
            end
        elseif i==7
            A{s}.avg = squeeze(mean(thisd(s,:,:,:),2));
        end
    end

    B = stai.trait;

    % Do t test
    cfg = [];
    cfg.neighbours = neighbours;
    cfg.channel ={'EEG','-M1','-M2'};
%     cfg.channel = {'C1','C2','Cz','CPz','CP1','CP2','Pz','P1','P2'};

    cfg.method = 'montecarlo';
    cfg.tail = 0;
    cfg.alpha = .025;

    cfg.correctm = 'cluster';
    cfg.clusteralpha = .05;
    cfg.minnbchan = 2;
    cfg.numrandomization = 500;

    cfg.statistic = 'ft_statfun_correlationT'; % 'indepsamplesregrT'
    cfg.design = B;
    cfg.ivar = 1;

    cfg.latency = [0 3];

    cfg.avgoverchan = 'no';
    cfg.avgovertime = 'no';

    stat = ft_timelockstatistics(cfg, A{:});

    figure
    imagesc(1-stat.prob);
    colormap('hot')
    xlabel('Time')
    ylabel('Channels')
    set(gca,'ytick',1:length(stat.label));
    set(gca,'yticklabels',stat.label);
    title(['Lowest p = ' num2str(min(stat.prob(:)))])
    set(gcf,'position',[440 -141 484 957])
    drawnow;

end

%% Onsets

mindur = 10; % how many samples of significant difference there needs to be for it to be an "onset"

onsetmethod = 'deviation'; % 'ttest': use t test to find pvalues < 0.05, 'deviation': use z scoring method
chanselection = {}; % {} to ignore and do all

onsets = [];
for s = 1:N

    % Load data
    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end

    [SL,RL,T,idx] = select_data(fullfile(dir_data,'newepoched',[schar '_FT_epoched_bc' num2str(baselineCorrect(b)) '.mat']),neighbours,standardType{st},datatype);
    
    if ~isempty(chanselection)
        cfg = [];
        cfg.channel = chanselection;
        cfg.avgoverchan = 'yes';
        SL = ft_selectdata(cfg,SL);
    end

    fidx = find(idx);
    nTrls = length(fidx);
    nChan = length(SL.label);

    % Get data matrix & rate of change
    X = SL.time{1};
    xidx = X<=3;

    roc = [];
    d = [];
    for trl = 1:nTrls
        if trl==1
            roc = nan(nTrls,nChan,sum(xidx));
        end
        for chan = 1:nChan
            tmp = SL.trial{fidx(trl)}(chan,xidx);
            tmp = movmean(tmp,10); % smooth
            d(trl,chan,:) = tmp;
            roc(trl,chan,2:end) = diff(tmp);
        end
    end

    roc = d;

    % Find first time that the rate of change differs between deviants vs standards
    for chan = 1:nChan

        switch onsetmethod
            case 'ttest'
                pvals = [];
                for t = 2:size(roc,3)
        
                    % neutral PE
                    y1 = abs(squeeze(roc(T.Condition(idx)==1,chan,t)));
                    y2 = abs(squeeze(roc(T.Condition(idx)==2,chan,t)));
                    [~,p] = ttest2(y1',y2','tail','left');
                    pvals(1,t-1) = p;
        
                    % fearful PE
                    y1 = abs(squeeze(roc(T.Condition(idx)==3,chan,t)));
                    y2 = abs(squeeze(roc(T.Condition(idx)==4,chan,t)));
                    [h,p] = ttest2(y1',y2','tail','left');
                    pvals(2,t-1) = p;
                end
                pvals = movmean(pvals,10);
        
                pvals = double(pvals < .05);
                theseonsets = pvals(:,1);
                for t = 2:size(pvals,2)
                    theseonsets(:,t) = theseonsets(:,t-1) + pvals(:,t);
                    theseonsets(pvals(:,t)==0,t) = 0;
                end
            case 'deviation'

                diffwave = [];

                y1 = squeeze(mean(roc(T.Condition(idx)==1,chan,2:end)));
                y2 = squeeze(mean(roc(T.Condition(idx)==2,chan,2:end)));
                diffwave(1,:) = zscore(y2-y1);

                y1 = squeeze(mean(roc(T.Condition(idx)==3,chan,2:end)));
                y2 = squeeze(mean(roc(T.Condition(idx)==4,chan,2:end)));
                diffwave(2,:) = zscore(y2-y1);

                diffwave = double(abs(diffwave) > 1);

                theseonsets = diffwave(:,1);
                for t = 2:size(diffwave,2)
                    theseonsets(:,t) = theseonsets(:,t-1) + diffwave(:,t);
                    theseonsets(diffwave(:,t)==0,t) = 0;
                end

        end

        for c = 1:2
            tmp = find(theseonsets(c,:)==mindur,1,'first');
            if ~isempty(tmp)
                onsets(s,chan,c) = tmp;
            else
                onsets(s,chan,c) = NaN;
            end
        end
    end
end

% Plot
figure
for chan = 1:nChan

    if isempty(chanselection)
        subplot(8,8,chan)
    end

    thisy = squeeze(onsets(:,chan,:));

    if sum(all(~isnan(thisy),2))>2
        m = nanmean(thisy);
        sem = nanstd(thisy)/sqrt(size(thisy,1));
        
        for c = 1:2
            plot(repmat(c,2,1),[m(c)+sem(c) m(c)-sem(c)],'k'); hold on
        end
        bar(m)
    
        [~,p] = ttest(thisy(:,1),thisy(:,2),'tail','right');
    
        title([SL.label{chan} ' p = ' num2str(round(p,3)) ', n = ' num2str(sum(~isnan(thisy)))])
    end
end

if isempty(chanselection)

    figure
    thisy = squeeze(mean(onsets(:,ismember(SL.label,{'CPz','Cz','FCz','Pz'}),:),2));
    m = nanmean(thisy);
    sem = nanstd(thisy)/sqrt(size(thisy,1));
    for c = 1:2
        plot(repmat(c,2,1),[m(c)+sem(c) m(c)-sem(c)],'k'); hold on
    end
    bar(m)
    [~,p] = ttest(thisy(:,1),thisy(:,2),'tail','right');
    title([SL.label{chan} ' p = ' num2str(round(p,3)) ', n = ' num2str(sum(~isnan(thisy)))])

    figure
    h = histfit(thisy(:,1));
    p1 = h(2).YData;
    x1 = h(2).XData;
    h = histfit(thisy(:,2));
    p2 = h(2).YData;
    x2 = h(2).XData;
    clf
    plot(x1,p1); hold on
    plot(x2,p2); hold on

end