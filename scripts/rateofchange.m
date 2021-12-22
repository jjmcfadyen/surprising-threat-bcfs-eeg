% Rate of change

onsetsamplemin = 2; % minimum sample time for an "onset" event
onsetsamplegap = 10; % time needed between onset events (in samples)
driftwindow = [0.5 0.1];

Y = [];
ROC = [];
onsets = [];
rlonsets = [];
trialdrifts = [];
alltrialonsets = cell(N,4);
allrltrialonsets = cell(N,4);
alltrialdrifts = cell(N,4);
rts = nan(N,4);
missingdata = nan(N,64);
onsetdata = [];
rlonsetdata = [];
rtonsetdata = [];
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

    % See when each trial & each channel exceeds a certain threshold
    trlonsets = nan(nTrls,nChan);
    rltrlonsets = nan(nTrls,nChan);
    parfor trl = 1:nTrls
        thisrt = T.RT(fidx(trl));
        thisx = X(xidx);
        thisxidx = thisx>0 & thisx<3; % thisxidx = thisx>0 & thisx<(thisrt-driftwindow(1));
        thisx = thisx(thisxidx);
        disp(['Trial ' num2str(trl)])
        for chan = 1:nChan

%             thisdata = squeeze(d(trl,chan,:));
%             thisdata = abs(zscore(thisdata(thisxidx)));
% 
%             thisroc = squeeze(roc(trl,chan,:));
%             thisroc = movmean(thisroc,SL.fsample*0.05); % smooth with 50 ms moving average (only if ROC)
%             thisroc = abs(zscore(thisroc(thisxidx)));

%             thisonset = find(cumsum(thisdata>1.28)>=onsetsamplemin,1,'first');
            
%             if ~isempty(thisonset)
%                 trlonsets(trl,chan) = thisx(thisonset-onsetsamplemin+1);
%             end

            thisonset = [];
            thisresid = [];
            tp = thisx(findMin(0,thisx)+1);
            tplist = [tp];
            while round(tp,3) <= thisx(end)
                thisdata = squeeze(d(trl,chan,thisx<=tp));
                [o,r] = findchangepts(thisdata,'statistic','mean');
                if isempty(o)
                    o = NaN;
                    r = NaN;
                end
                thisonset = [thisonset o];
                thisresid = [thisresid r];
                try
                    tp = thisx(findMin(tp,thisx)+10);
                    tplist = [tplist tp];
                catch
                    tp = Inf;
                end
            end
            
%             thisonset = mode(thisonset(tplist<(thisrt-driftwindow(1)))); % find the most common onset point prior to response time (plus motor preparation window)
            thisonset = mode(thisonset(tplist<thisrt));
%             thisonset = findchangepts(squeeze(d(trl,chan,thisx>0 & thisx<thisrt)),'statistic','mean'); % find most "significant" change point in window from zero to response onset

            if ~isempty(thisonset)
                trlonsets(trl,chan) = thisx(thisonset(1));
                rltrlonsets(trl,chan) = -(thisrt - thisx(thisonset));
            end
        end
    end
    disp([':::::::::::: ' num2str(round(mean(mean(isnan(trlonsets)))*100,2)) '% of trials with no identifiable onset time'])
    disp([':::::::::::: Mean onset time = ' num2str(round(nanmean(trlonsets(:)),3)) ' seconds'])
    disp([':::::::::::: Mean onset time = ' num2str(round(nanmean(rltrlonsets(:)),3)) ' seconds response locked)'])
    disp([':::::::::::: Mean no. of onsets after response = ' num2str(nanmean(rltrlonsets(:)>=0))])

    missingdata(s,:) = mean(isnan(trlonsets));

    for c = 1:4
        onsets(s,c,:) = nanmean(trlonsets(T.Condition(idx)==c,:));
        rlonsets(s,c,:) = nanmean(rltrlonsets(T.Condition(idx)==c,:));
        alltrialonsets{s,c} = trlonsets(T.Condition(idx)==c,:);
        allrltrialonsets{s,c} = rltrlonsets(T.Condition(idx)==c,:);
    end

    % Gather data to plot, organised by onset time
    thisx = X(xidx);
    for chan = 1:nChan
        q = quantile(trlonsets(:,chan),4);
        for i = 1:3

            % get index of trials with this onset time (fast, medium, or slow)
            qidx = trlonsets(:,chan)>q(i) & trlonsets(:,chan)<q(i+1);
            fqidx = find(qidx);

            % select stimulus-locked data for this onset group and average it
            onsetdata(s,i,chan,:) = squeeze(mean(d(qidx,chan,:)));

            % log average onset times for this channel/quantile
            rtonsetdata(s,i,chan,1) = mean(trlonsets(qidx,chan));

        end
    end

    % Gather response-locked data to plot, organised by onset time
    thisx = X(xidx);
    theserts = T.RT(idx);
    for chan = 1:nChan
        q = fliplr(quantile(rltrlonsets(:,chan),4)); % fastest to slowest
        for i = 1:3

            % get index of trials with this onset time (fast, medium, or slow)
            qidx = rltrlonsets(:,chan)<q(i) & rltrlonsets(:,chan)>q(i+1);
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
            rlonsetdata(s,i,chan,:) = thisrlmean;

            % log average onset times for this channel/quantile
            rtonsetdata(s,i,chan,2) = mean(rltrlonsets(qidx,chan));

        end
    end

%     avcon = [];
%     for c = 1:4
%         avcon(c,:,:) = mean(d(T.Condition(idx)==c,:,:));
%     end
% 
%     thisx = X(xidx);
%     for c = 1:4
%         for chan = 1:nChan
% 
%             thisonset = [];
%             thisresid = [];
%             tp = thisx(findMin(0.5,thisx));
%             while tp < 3
%                 thisdata = squeeze(avcon(c,chan,thisx>0 & thisx<=tp));
%                 [o,r] = findchangepts(thisdata);
%                 thisonset = [thisonset o];
%                 thisresid = [thisresid r];
%                 tp = thisx(findMin(tp,thisx)+1);
%             end
%             
%             ri = findchangepts(thisresid);
%             thisonset = thisonset(ri);
% 
%             onsets(s,c,chan) = thisx(thisonset + findMin(0,thisx));
%         end
%     end

%     % See when all channels differ significantly from baseline
%     tmpx = X(xidx);
%     trlonsets = nan(1,nTrls);
%     for trl = 1:nTrls
%         thisbaseline = squeeze(nanmean(abs(roc(trl,:,tmpx<=0)),3));
%         tp = find(tmpx==0);
%         pvals = [];
%         while true
%             tp = tp+1;
%             if tp>=size(roc,3)
%                 break
%             end
%             [~,p] = ttest(abs(squeeze(roc(trl,:,tp))),thisbaseline,'tail','right');
%             pvals = [pvals p];
%             if length(pvals)>=onsetsamplemin
%                 if all(pvals(end-onsetsamplemin+1:end) < .05)
%                     trlonsets(trl) = tmpx(tp-onsetsamplemin);
%                     break
%                 end
%             end
%             if isnan(trlonsets(trl))
%                 trlonsets(trl) = tmpx(find(pvals==min(pvals),1,'first'));
%             end
%         end
%     end
% 
%     for c = 1:4
%         onsets(s,c) = nanmean(trlonsets(T.Condition(idx)==c));
%         alltrialonsets{s,c} = trlonsets(T.Condition(idx)==c);
%     end

%     % Do stats against zero OR baseline
%     nTime = size(d,3);
%     thisx = X(xidx);
%     thisx = thisx(thisx>0);
%     for i = 1:2
% 
%         theseonsets = nan(nChan,4,2);
%         for chan = 1:nChan
% 
%             if i==1
%                 thisnull = zeros(size(d,1),2);
%             elseif i==2
%                 thisnull = [squeeze(mean(d(:,chan,X<=0),3)) squeeze(nanmean(roc(:,chan,X<=0),3))];
%             end
% 
%             % Get pvalue of each condition compared to a null (zeros or baseline period), at each time point
%             pvals = nan(4,nTime,2); % condition, time point, data or rate of change
%             for c = 1:4
% 
%                 % get average waveform for this condition & this channel
%                 cidx = T.Condition(idx)==c;
%                 thisy = squeeze(d(cidx,chan,:));
%                 thisy(:,:,2) = squeeze(roc(cidx,chan,:));
% 
%                 % reapply baseline correction on the actual data (not the rate of change)
%                 thisy(:,:,1) = thisy(:,:,1) - mean(thisy(:,X(xidx)>= -0.05 & X(xidx)<=0,1),2);
% 
%                 for t = 1:nTime
%                     thisy = [squeeze(d(cidx,chan,t)) squeeze(roc(cidx,chan,t))];
%                     [~,pvals(c,t,1)] = ttest(thisy(:,1),thisnull(cidx,1));
%                     if t>1
%                         [~,pvals(c,t,2)] = ttest(thisy(:,2),thisnull(cidx,2));
%                     end
%                 end
%             end
%             pvals = pvals(:,X(xidx)>0,:); % don't look at pvalues in baseline
% 
%             % if the channel isn't dead (all zeros/nans)
%             if ~all(isnan(pvals(:)))
% 
%                 % Correct for multiple comparisons (FDR) across all timepoints (per channel/dataset)
%                 P = squeeze(pvals(:,:,1));
%                 Pd = reshape(mafdr(P(:)),size(P,1),size(P,2)) < .05;
% 
%                 P = squeeze(pvals(:,:,2)); 
%                 Pr = reshape(mafdr(P(:)),size(P,1),size(P,2)) < .05;
% 
%                 % See if there is a significant onset for each condition 
%                 % --- data against null
%                 if all(any(Pd,2))
% 
%                     % Move through time window and find onsets of significant clusters (minimum = onsetsamplemin)
%                     pcount = double(Pd);
%                     for t = 2:size(Pd,2)
%                         pcount(:,t) = pcount(:,t-1) + Pd(:,t);
%                         pcount(Pd(:,t)==0,t) = 0;
%                     end
% 
%                     for c = 1:4
%                         o = find(pcount(c,:)>=onsetsamplemin,1,'first');
%                         if ~isempty(o)
%                             theseonsets(chan,c,1) = thisx(o);
%                         end
%                     end
%                 end
%                 if all(any(Pr,2))
%                     
%                     % Move through time window and find onsets of significant clusters (minimum = onsetsamplemin)
%                     pcount = double(Pr);
%                     for t = 2:size(Pr,2)
%                         pcount(:,t) = pcount(:,t-1) + Pr(:,t);
%                         pcount(Pr(:,t)==0,t) = 0;
%                     end
% 
%                     for c = 1:4
%                         o = find(pcount(c,:)>=onsetsamplemin,1,'first');
%                         if ~isempty(o)
%                             theseonsets(chan,c,2) = thisx(o);
%                         end
%                     end
% 
%                 end
%             end
%         end
%         if i==1
%             onsets(s,:,:,:) = theseonsets;
%         elseif i==2
%             onsets_baseline(s,:,:,:) = theseonsets;
%         end
%     end

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
    'onsets','rlonsets','alltrialonsets','allrltrialonsets','rts','missingdata','onsetdata','rlonsetdata','rtonsetdata','onsetsamplemin','onsetsamplegap','driftwindow','labels','x');

%% Plot

chanselection = [];
chanselection{1} = {'PO4','PO8','O2','PO3','PO7','O1','Oz','POz'};
chanselection{2} = {'CPz','Cz','Pz'};

cmap = [243, 140, 255 ;
        178, 71, 255 ;
        66, 0, 213]/255;

% Data grouped by onset time
% figure
% for chan = 1:size(Y,3)
% 
%     subplot(8,8,chan)
% 
%     for i = 1:3
% 
%         y = squeeze(onsetdata(:,i,chan,:));
%         m = mean(y);
%         sem = std(y)/sqrt(size(y,1));
%         upper = m+sem;
%         lower = m-sem;
% 
%         patch([x fliplr(x)],[upper fliplr(lower)],cmap(i,:),'facealpha',.2,'edgealpha',0); hold on
%         plot(x,m,'color',cmap(i,:),'linewidth',1.3); hold on
%     end
%     
%     ax = gca;
%     for i = 1:3
%         plot(repmat(mean(rtonsetdata(:,i,chan,1)),2,1),ax.YLim,'color',cmap(i,:),'linewidth',1.2,'linestyle','--'); hold on
%     end
%     plot(x([1 end]),[0 0],'k:'); hold on
%     title(labels{chan})
%     xlim(x([1 end]))
% 
% end
% sgtitle('SL onset times (fast, medium, slow)')

figure
for chan = 1:length(chanselection) % specific channel groupings

    subplot(1,length(chanselection),chan)

    for i = 1:3

        y = squeeze(mean(onsetdata(:,i,ismember(labels,chanselection{chan}),:),3));
        m = mean(y);
        sem = std(y)/sqrt(size(y,1));
        upper = m+sem;
        lower = m-sem;

        patch([x fliplr(x)],[upper fliplr(lower)],cmap(i,:),'facealpha',.2,'edgealpha',0); hold on
        plot(x,m,'color',cmap(i,:),'linewidth',1.3); hold on
    end
    
    ax = gca;
    for i = 1:3
        plot(repmat(mean(rtonsetdata(:,i,chan,1)),2,1),ax.YLim,'color',cmap(i,:),'linewidth',1.2,'linestyle','--'); hold on
    end
    plot(x([1 end]),[0 0],'k:'); hold on
    title(strjoin(chanselection{chan}))
    xlim(x([1 end]))

end
sgtitle('SL onset times (fast, medium, slow)')

% -- response-locked
thisx = x - 2.9;

% figure
% for chan = 1:size(Y,3)
% 
%     subplot(8,8,chan)
% 
%     for i = 1:3
% 
%         y = squeeze(rlonsetdata(:,i,chan,:));
%         m = nanmean(y);
%         sem = nanstd(y)/sqrt(size(y,1));
%         upper = m+sem;
%         lower = m-sem;
% 
%         patch([thisx fliplr(thisx)],[upper fliplr(lower)],cmap(i,:),'facealpha',.2,'edgealpha',0); hold on
%         plot(thisx,m,'color',cmap(i,:),'linewidth',1.3); hold on
%     end
%     
%     ax = gca;
%     for i = 1:3
%         plot(repmat(nanmean(rtonsetdata(:,i,chan,2)),2,1),ax.YLim,'color',cmap(i,:),'linewidth',1.2,'linestyle','--'); hold on
%     end
%     tmpx = ~isnan(squeeze(nanmean(nanmean(rlonsetdata(:,:,chan,:)),2)));
%     tmpx = thisx([find(tmpx,1,'first') find(tmpx,1,'last')]);
%     plot(tmpx,[0 0],'k:'); hold on
%     plot([0 0],ax.YLim,'k:'); hold on
%     title(SL.label{chan})
%     xlim([-1.5 0.1])
% 
% end
% sgtitle('RL onset times (fast, medium, slow)')

figure
for chan = 1:length(chanselection) % specific channel groupings

    subplot(1,length(chanselection),chan)

    for i = 1:3

        y = squeeze(nanmean(rlonsetdata(:,i,ismember(labels,chanselection{chan}),:),3));
        m = nanmean(y);
        sem = nanstd(y)/sqrt(size(y,1));
        upper = m+sem;
        lower = m-sem;

        nanidx = isnan(upper) | isnan(lower);
        patch([thisx(~nanidx) fliplr(thisx(~nanidx))],[upper(~nanidx) fliplr(lower(~nanidx))],cmap(i,:),'facealpha',.2,'edgealpha',0); hold on
        plot(thisx,m,'color',cmap(i,:),'linewidth',1.3); hold on
    end
    
    ax = gca;
    for i = 1:3
        plot(repmat(nanmean(rtonsetdata(:,i,chan,2)),2,1),ax.YLim,'color',cmap(i,:),'linewidth',1.2,'linestyle','--'); hold on
    end
    tmpx = ~isnan(squeeze(nanmean(nanmean(rlonsetdata(:,:,chan,:)),2)));
    tmpx = thisx([find(tmpx,1,'first') find(tmpx,1,'last')]);
    plot(tmpx,[0 0],'k:'); hold on
    plot([0 0],ax.YLim,'k:'); hold on
    title(strjoin(chanselection{chan}))
    xlim([-1.5 0.1])

end
sgtitle('RL onset times (fast, medium, slow)')

% % Data itself
% 
% cmap = [0 232 255;
%         51 110 255;
%         255 186 48;
%         255 0 0]/255;
% 
% figure
% for chan = 1:size(Y,3)
% 
%     subplot(8,8,chan)
%     for c = 1:4
% 
%         y = squeeze(Y(:,c,chan,:));
%         % z score each subject?
%         m = mean(y);
%         sem = std(y)/sqrt(size(y,1));
%         upper = m+sem;
%         lower = m-sem;
% 
%         patch([x fliplr(x)],[upper fliplr(lower)],cmap(c,:),'facealpha',.2,'edgealpha',0); hold on
%         plot(x,m,'color',cmap(c,:),'linewidth',1.3); hold on
%     end
%     plot(x([1 end]),[0 0],'k:'); hold on
%     title(SL.label{chan})
%     xlim(x([1 end]))
%     
% end
% sgtitle('Data')
% 
% figure
% for chan = 1:size(ROC,3)
% 
%     subplot(8,8,chan)
%     for c = 1:4
% 
%         y = squeeze(ROC(:,c,chan,2:end));
%         m = mean(y);
%         sem = std(y)/sqrt(size(y,1));
%         upper = m+sem;
%         lower = m-sem;
% 
%         patch([x(2:end) fliplr(x(2:end))],[upper fliplr(lower)],cmap(c,:),'facealpha',.2,'edgealpha',0); hold on
%         plot(x(2:end),m,'color',cmap(c,:),'linewidth',1.3); hold on
%     end
%     plot(x([1 end]),[0 0],'k:'); hold on
%     title(SL.label{chan})
%     xlim(x([1 end]))
% end
% sgtitle('Rate of change')

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

    for c = 1:4

        nTrls = size(alltrialonsets{s,c},1);

        thistable = array2table([repmat(s,nTrls,1) repmat(c,nTrls,1)],'variablenames',{'Subject','Condition'});
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

        for chan = 1:length(chanselection)
            thistable.([changroupnames{1} 'onsetSL']) = mean(alltrialonsets{s,c}(:,ismember(labels,chanselection{chan})),2);
            thistable.([changroupnames{1} 'onsetRL']) = mean(allrltrialonsets{s,c}(:,ismember(labels,chanselection{chan})),2);
            thistable.([changroupnames{1} 'drift']) = mean(alltrialdrifts{s,c}(:,ismember(labels,chanselection{chan})),2);
        end

        thistable.RT = T.RT(idx & T.Condition==c);

        thistable.nondecision = nan(nTrls,1);
        thistable.drift = nan(nTrls,1);
        thistable.boundary = nan(nTrls,1);

        if any(ismember(ddm.Subject,subjects(s)))
            thistable.nondecision = repmat(table2array(ddm(ismember(ddm.Subject,subjects(s)),ismember(ddm.Properties.VariableNames,['nondecision_C' num2str(c)]))),nTrls,1);
            thistable.drift = repmat(table2array(ddm(ismember(ddm.Subject,subjects(s)),ismember(ddm.Properties.VariableNames,['drift_C' num2str(c)]))),nTrls,1);
            thistable.boundary = repmat(table2array(ddm(ismember(ddm.Subject,subjects(s)),ismember(ddm.Properties.VariableNames,['boundary_C' num2str(c)]))),nTrls,1);
        end
    end

    thistable.Anxiety = repmat(stai.trait(s),size(thistable,1),1);

    % add to table
    L = [L; thistable];

end

writetable(L,'D:\bCFS_EEG_Reanalysis\results\erponsets_alltrials.csv');


%%
% Onsets & drifts
ddm = readtable('D:\bCFS_EEG_Reanalysis\results\ddmparams_nozscore.csv'); 
ddm.Subject(ddm.Experiment==2) = ddm.Subject(ddm.Experiment==2)-31; % 31 subjects in experiment 1
ddm = ddm(ddm.Experiment==2,:);
ddm = ddm(ismember(ddm.Subject,subjects),:);

% ter = [ddm.nondecision_EN ddm.nondecision_UN ddm.nondecision_EF ddm.nondecision_UF];
% v = [ddm.drift_EN ddm.drift_UN ddm.drift_EF ddm.drift_UF];

% cmap = [141, 34, 255
%     229, 131, 255]/255;
% for i = 1:2
%     
%     figure
% 
%     if i==1
%         thisddm = ter;
%     elseif i==2
%         thisddm = v;
%     end
%     thisddm = thisddm(:);
% 
%     for chan = 1:size(Y,3)
%     
%         thiseeg = squeeze(Y(ismember(subjects,ddm.Subject),:,chan,:));
%         thiseeg = reshape(thiseeg,size(thiseeg,1)*size(thiseeg,2),size(thiseeg,3));
% 
%         subplot(8,8,chan)
%         for c = 1:2 % low, high
%     
%             if c==1
%                 cidx = thisddm<median(thisddm);
%             elseif c==2
%                 cidx = thisddm>median(thisddm);
%             end
% 
%             y =  thiseeg(cidx,:);
%             m = mean(y);
%             sem = std(y)/sqrt(size(y,1));
%             upper = m+sem;
%             lower = m-sem;
%     
%             patch([x fliplr(x)],[upper fliplr(lower)],cmap(c,:),'facealpha',.2,'edgealpha',0); hold on
%             plot(x,m,'color',cmap(c,:),'linewidth',1.3); hold on
%         end
%         plot(x([1 end]),[0 0],'k:'); hold on
%         title(SL.label{chan})
%         
%     end
%     if i==1
%         sgtitle('Nondecision time')
%     elseif i==2
%         sgtitle('Drift rate')
%     end
% end
% 
% zonsets = onsets;
% for s = 1:N
%     for chan = 1:nChan
%         zonsets(s,:,chan) = zscore(squeeze(onsets(s,:,chan)));
%     end
% end
% 
% % remove onset outliers
% ronsets = onsets;
% rzonsets = zonsets;
% for s = 1:N
%     for chan = 1:nChan
%         tmp = [];
%         for c = 1:4
%             tmp = [tmp; alltrialonsets{s,c}(:,chan) ones(size(alltrialonsets{s,c},1),1)*c];
%         end
%         ridx = (tmp(:,1) - nanmean(tmp(:,1))) / nanstd(tmp(:,1));
%         ridx = double(abs(ridx) > 3);
%         ridx(isnan(ridx)) = 1;
%         for c = 1:4
%             ronsets(s,c,chan) = nanmean(tmp(tmp(:,2)==c & ~ridx,1));
%         end
%         rzonsets(s,:,chan) = zscore(squeeze(ronsets(s,:,chan)));
%     end
% end

cmap = [0 232 255;
        51 110 255;
        255 186 48;
        255 0 0]/255;

% figure
% for chan = 1:64
%     subplot(8,8,chan)
%     thisy = squeeze(onsets(:,:,chan));
%     thisy = zscore(thisy')';
%     m = mean(thisy);
%     sem = std(thisy)/sqrt(size(thisy,1));
%     upper = m+sem;
%     lower = m-sem;
%     for c = 1:4
%         plot([c c],[lower(c) upper(c)],'k','linewidth',1.3); hold on
%         scatter(c,m(c),70,'markerfacecolor',cmap(c,:),'markeredgecolor','k'); hold on
%     end
%     [~,p] = ttest(thisy(:,2)-thisy(:,1),thisy(:,4)-thisy(:,3));
%     title([SL.label{chan} ', p = ' num2str(round(p,3))])
% end
% sgtitle('ROC onsets')

for j = 1%:2
    for i = 1:2
    
        if i==1
            chanselection = {'PO4','PO8','O2','PO3','PO7','O1','Oz','POz'};
        elseif i==2
            chanselection = {'CPz','Cz','Pz'};
        end
    
        if j==1 % stimulus-locked onsets
            thisy = squeeze(mean(onsets(:,:,ismember(labels,chanselection)),3));
        elseif j==2 % response-locked onsets
            thisy = squeeze(mean(rlonsets(:,:,ismember(labels,chanselection)),3));
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
        title([strjoin(chanselection,' ') ', p = ' num2str(round(p,3))])
        xlim([0 5])
        set(gca,'ticklength',[0 0])
        set(gcf,'position',[715 332 296 336])
    
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
        for j = 1:4
    
            figure
    
            if j==1
                thisparam = 'nondecision';
            elseif j==2
                thisparam = 'drift';
            elseif j==3
                thisparam = 'boundary';
            elseif j==4
                thisparam = 'rt';
            end
    
            if j<=3
                thiseeg = squeeze(mean(onsets(ismember(subjects,ddm.Subject),:,ismember(SL.label,chanselection)),3));
                thisddm = table2array(ddm(:,contains(ddm.Properties.VariableNames,thisparam)));
            elseif j==4
                thiseeg = squeeze(mean(onsets(:,:,ismember(SL.label,chanselection)),3));
                thisddm = rts;
                
            end
            scatter(thiseeg(:),thisddm(:),'markerfacecolor','k','markeredgecolor','none','markerfacealpha',.5); hold on
            coeff = polyfit(thiseeg(:),thisddm(:),1);
            thisline = polyval(coeff,thiseeg(:));
            plot(thiseeg(:),thisline,'k'); hold on
            [r,p] = corr(thiseeg(:),thisddm(:));
            title([thisparam ': r=' num2str(round(r,3)) ', p=' num2str(round(p,3))])
        end
    end
end


% thisy = [];
% for s = 1:N
%     suby = [];
%     for c = 1:4
%         suby = [suby; alltrialonsets{s,c}' ones(length(alltrialonsets{s,c}),1)*c];
%     end
%     suby(:,1) = zscore(suby(:,1));
% %     suby(abs(suby(:,1))>3,:) = []; % remove outliers
%     for c = 1:4
%         thisy(s,c) = mean(suby(suby(:,2)==c,1));
%     end
% end
% 
% figure
% for c = 1:4
%     binwidth = diff(linspace(min(thisy(c,:)),max(thisy(c,:)),20));
%     binwidth = binwidth(1);
% 
%     [x,y] = beeswarm(thisy(:,c),binwidth,0.15);
% 
%     q = quantile(y,[.25 .75]);
%     patch([c-.15 c-.15 c+.15 c+.15],[q(1) q(2) q(2) q(1)],cmap(c,:),'facealpha',.25,'edgecolor',cmap(c,:)); hold on
%     plot([c c],[q(2) max(thisy(:,c))],'color',cmap(c,:),'linewidth',1.2); hold on
%     plot([c c],[q(1) min(thisy(:,c))],'color',cmap(c,:),'linewidth',1.2); hold on
% 
%     scatter(x+c,y,30,'markerfacecolor',cmap(c,:),'markeredgecolor','none','markerfacealpha',.5); hold on
% 
%     m = mean(y);
%     sem = std(y)/sqrt(size(y,1));
%     upper = m+sem;
%     lower = m-sem;
% 
%     plot([c c],[lower upper],'k','linewidth',1.3); hold on
%     scatter(c,m,70,'markerfacecolor',cmap(c,:),'markeredgecolor','k'); hold on
% 
% end
% xlim([0 5])
% set(gca,'ticklength',[0 0])
% sgtitle('ROC onsets')
% 
% 
% tmpx = X(xidx);
% avonsets = nan(N,4);
% for s = 1:N
%     for c = 1:4
%         thisbaseline = abs(squeeze(nanmean(Y(s,c,:,tmpx<=0),4)));
%         tp = find(tmpx==0);
%         pvals = [];
%         while true
%             tp = tp+1;
%             if tp>=size(Y,3)
%                 break
%             end
%             [~,p] = ttest(abs(squeeze(Y(s,c,:,tp))),thisbaseline,'tail','right');
%             pvals = [pvals p];
%             if length(pvals)>=onsetsamplemin
%                 if all(pvals(end-onsetsamplemin+1:end) < .05)
%                     avonsets(s,c) = tmpx(tp-onsetsamplemin);
%                     break
%                 end
%             end
%             if isnan(avonsets(s,c))
%                 avonsets(s,c) = tmpx(find(pvals==min(pvals),1,'first'));
%             end
%         end
%     end
% end
% avonsets = zscore(avonsets')';




% 
% 
% figure
% for chan = 1:64
%     subplot(8,8,chan)
%     thisy = squeeze(abs(trialdrifts(:,:,chan)));
%     m = mean(thisy);
%     sem = std(thisy)/sqrt(size(thisy,1));
%     upper = m+sem;
%     lower = m-sem;
%     for c = 1:4
%         plot([c c],[lower(c) upper(c)],'k','linewidth',1.3); hold on
%         scatter(c,m(c),70,'markerfacecolor',cmap(c,:),'markeredgecolor','k'); hold on
%     end
%     title([SL.label{chan}])
% end
% sgtitle('Drift -800 to -200 ms pre-response onset')
% 
% 
% for i = 1:2
%     figure
%     for chan = 1:64
%     
%         subplot(8,8,chan)
%     
%         thisy = squeeze(onsets(ismember(subjects,ddm.Subject),chan,:,1));
%         thisy = zscore(thisy')';
% 
%         if i==1
%             thisx = ter;
%         elseif i==2
%             thisx = v;
%         end
%     
%         thisy = thisy(:);
%         thisx = thisx(:);
%     
%         nanidx = isnan(thisy);
%         thisy = thisy(~nanidx);
%         thisx = thisx(~nanidx);
%     
% %         [rho,p] = corr(thisx,thisy);
% %         scatter(thisy,thisx); hold on
% %         coeff = polyfit(thisy,thisx,1);
% %         fitline = polyval(coeff,thisy);
% %         plot(thisy,fitline,'k');
%         
%         m = mean([thisy(thisx<median(thisx)) thisy(thisx>median(thisx))]);
%         sd = std([thisy(thisx<median(thisx)) thisy(thisx>median(thisx))]);
%         
%         bar(m); hold on
%         for j = 1:2
%             plot([j j],[m(j)-(sd(j)/2) m(j)+(sd(j)/2)],'k');
%         end
%         [h,p] = ttest(thisy(thisx<median(thisx)),thisy(thisx>median(thisx)));
% 
%         title([SL.label{chan} ': p = ' num2str(round(p,3))])
%     end
%     if i==1
%         sgtitle('Onset vs Nondecision time')
%     elseif i==2
%         sgtitle('Onset vs Drift rate')
%     end
% end

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
