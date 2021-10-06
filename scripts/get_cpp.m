function get_cpp(subjects,st,baselineType)
% st = standard types (1 = all, 2 = first, 3 = last)
% baselineType = 'stim' or 'resp'

dir_data = 'D:\bCFS_EEG_Reanalysis\data\Exp2\eeg\preprocessed\unfiltered';
dir_save = 'D:\bCFS_EEG_Reanalysis\results\cpp';

N = length(subjects);

%% Select CPP channel for each participant

datatype = 'SL'; % 'SL' or 'RL'
rtlims = [-1 0]; % how long for the epochs to be (at least) before the response
datamethod = 'GLM'; % 'GLM' (to do a glm between ALL channels and RT),
                    % 'single' (to get single central-parietal channel with highest correlation between positive peak and RT)
                    % 'average' (to get an average of central-parietal channels with positive peaks at RT)

X = [];
Y = [];
xyinfo = []; % col 1 = subject, col 2 = trial, to describe the first dimension of X and Y (pooled trials across participants)
if ~strcmp(datamethod,'GLM')
    channelselection = array2table([subjects' nan(N,2)],'variablenames',{'Subject','Channel','R'});
    channelselection.Channel = cell(N,1);
    collated = [];
    collated.eeg = [];
    collated.rt = [];
else
    rtbetas = [];
    varexplained = [];
    betaplot = [];
end
for s = 1:N

    subject = subjects(s);
    if subject < 10
        schar = ['S0' num2str(subject)];
    else
        schar = ['S' num2str(subject)];
    end
    
    disp(['Reading in response-locked data for ' schar '...'])
    
    % Get data
    switch baselineType
        case 'stim'
            thisbc = '1';
        case 'resp'
            thisbc = '0';
    end
    [SL,RL,T,idx] = clean_data(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' st 'Standards_bc' thisbc '.mat']));
    
    switch datatype
        case 'SL'
            data = SL;
        case 'RL'
            data = RL;
    end
    
    fidx = find(idx);
    nTrls = sum(idx);
    
%     % Downsample
%     cfg = [];
%     cfg.resamplefs = 100;
%     data = ft_resampledata(cfg,data);
    
    % Select trials & channels & time-windows of interest
    cfg = [];
    cfg.trials = fidx;
    cfg.channel = {'EEG','-M1','-M2'};
    thisdata = ft_selectdata(cfg,data);
    
    chanidx = ismember(data.label,thisdata.label);
    nChan = sum(chanidx);
    
    if strcmp(datatype,'RL') && strcmp(datamethod,'GLM')
        
        tlengths = nan(nTrls,2);
        for trl = 1:nTrls
            tlengths(trl,:) = thisdata.time{trl}([1 end]);
        end
        
        % exclude any trials that don't have at least 0 second post-response time or that have less than 1 second before response
        excludeIdx = tlengths(:,2) < rtlims(2) | tlengths(:,1) > rtlims(1);
        if any(excludeIdx)
            disp(['!!! Removing ' num2str(sum(excludeIdx)) ' trials (' num2str(round(mean(excludeIdx)*100,2)) '%) due to trial length'])
            
            idx(fidx(excludeIdx)) = 0;
            fidx(excludeIdx,:) = [];
            nTrls = sum(idx);
            
            cfg = [];
            cfg.trials = find(~excludeIdx);
            thisdata = ft_selectdata(cfg,thisdata);            
        end
        
        tlengths = nan(nTrls,2);
        for trl = 1:nTrls
            tlengths(trl,:) = thisdata.time{trl}([1 end]);
        end
        
        % crop the data
        cfg = [];
        cfg.toilim = [max(tlengths(:,1)) min(tlengths(:,2))];
        thisdata = ft_redefinetrial(cfg,thisdata);
        
    end
    
    avdata = ft_timelockanalysis([],thisdata);
    nTime = length(avdata.time);
    
    % GLM
    switch datamethod
        case 'GLM'
        
            dm = zeros(nTime*nTrls,nTrls+1);
            dm(:,1) = reshape(repmat(T.RT(idx)',nTime,1),nTrls*nTime,1); % RT
%             dm(:,1) = reshape(repmat(T.Condition(idx)>2',nTime,1),nTrls*nTime,1); % Emotion (fear=1, neutral=0)
            tmp = repmat({ones(nTime,1)},1,nTrls);
            dm(:,2:end) = blkdiag(tmp{:});

            thisy = nan(nTime*nTrls,nChan);
            cc = 1;
            for trl = 1:nTrls
                thisidx = cc:(cc+nTime-1);
                thisy(thisidx,:) = thisdata.trial{trl}';
                cc = thisidx(end)+1;
            end

            B = pinv(dm)*thisy;
            E = (thisy - dm*B);
            vary = sum(var(dm*B));
            vare = sum(var(E));
            varexplained(s,1) = vary / (vary + vare); % from here: https://web.stanford.edu/group/vista/cgi-bin/wiki/index.php/GLM
            
            rtbetas(s,:) = B(1,:);

            thisplot = avdata;
            thisplot.time = 0;
            thisplot.avg = B(1,:)';
            thisplot = rmfield(thisplot,{'var','sem','dof'});
            
            betaplot{s} = thisplot;

            figure

            subplot(1,2,1)
            cfg = [];
            cfg.layout = 'biosemi64.lay';
%             cfg.xlim = 0;
            ft_topoplotER(cfg,thisplot);
            colormap(colours(100,'viridis'))

            newy = [];
            cc = 1;
            for trl = 1:nTrls
                thisidx = cc:(cc+nTime-1);
                newy(trl,:) = thisy(thisidx,:)*B(1,:)';
                cc = thisidx(end)+1;
            end

            cmap = [0 232 255;
                51 110 255;
                255 186 48;
                255 0 0]/255;

            subplot(1,2,2)
            for c = 1:4
                plot(avdata.time,mean(newy(T.Condition(idx)==c,:)),'linewidth',1.3,'color',cmap(c,:)); hold on
            end
            sgtitle(schar)
            drawnow
            
            % Add to group variables
            switch datatype
                case 'SL'
                    X = [X; avdata.time];
                    if isempty(Y)
                        Y = newy;
                    else
                        Y = [Y; newy];
                    end
                case 'RL'
                    X{s} = avdata.time;
                    Y{s} = newy;
            end
    
        % See which mode covaries the most with RT
        case 'single'

            xlength = length(SL.time{1}); % use first trial to see how many samples are in all trials
            thismin = min([min(T.RT(fidx)) 1]); % only look 1 second before RT (or, if there are RTs shorter than this, then take the shortest RT)

            % smooth each channel's data using loess filter using 30% of length
            labelidx = ~contains(SL.label,{'M1','M2','Nz','SNz','LH','RH','LV','UV','Status'});
            nChan = sum(labelidx);
            
            allsmooth = nan(nTrls,nChan,xlength);
            for trl = 1:nTrls
                allsmooth(trl,:,:) = SL.trial{fidx(trl)}(labelidx,:);
            end

            disp('Smoothing...')
            for chan = 1:size(allsmooth,2)
                allsmooth(:,chan,:) = smoothdata(squeeze(allsmooth(:,chan,:)),2,'loess',round(0.3*xlength));
            end

            disp('Finding peaks...')
            pklocs = nan(nTrls,nChan);
            for trl = 1:nTrls

                thisrt = T.RT(fidx(trl));
                origx = SL.time{fidx(trl)};
                thisy = SL.trial{fidx(trl)}(labelidx,:);

                xidx = origx >= (thisrt-thismin) & origx <= thisrt+.2;
                thisx = origx(xidx);

                if size(thisy,1) ~= nChan
                    error([num2str(size(thisy,1)) ' channels detected instead of 64'])
                end

                smoothed = squeeze(allsmooth(trl,:,:));

                tmppklocs = nan(1,nChan);
                for chan = 1:nChan
                    [amp,loc] = findpeaks(smoothed(chan,xidx));
                    if isempty(amp) || smoothed(chan,findMin(origx,thisrt)) < 0
                        tmppklocs(1,chan) = NaN; % no positive peak, or amplitude is negative at time of response
                    else
                        tmppklocs(1,chan) = thisx(loc(amp==max(amp)));
                    end
                end
                pklocs(trl,:) = tmppklocs;

            end

            theselabels = SL.label(labelidx);
            chanidx = find(ismember(theselabels,{'C1','C2','Cz','CP1','CP2','CPz','P1','P2','Pz'}));
            rtcorr = nan(length(chanidx),1);
            rtdist = nan(length(chanidx),1);
            for chan = 1:length(chanidx)
                thisx = pklocs(:,chanidx(chan));
                thisx(isnan(thisx)) = 0;
                rtcorr(chan,1) = corr(thisx,T.RT(idx),'rows','pairwise');
                rtdist(chan,1) = mean((T.RT(idx) - thisx) .^ 2);
            end
            bestCombo = pdist2([rtcorr rtdist],[max(rtcorr) min(rtdist)]);
            bestChan = chanidx(bestCombo==min(bestCombo));

            channelselection.R(channelselection.Subject==subjects(s)) = rtcorr(1);
            channelselection.Channel{channelselection.Subject==subjects(s)} = theselabels{bestChan};

            % concatenate subject data for plotting RT correlation later
            [sortedrt,sortidx] = sort(T.RT(idx));

            sortedeeg = [];
            for trl = 1:nTrls
                sortedeeg(trl,:) = SL.trial{fidx(trl)}(bestChan,:); 
            end
            sortedeeg = sortedeeg(sortidx,:);

            collated.eeg = [collated.eeg; sortedeeg];
            collated.rt = [collated.rt; sortedrt];

            figure
            imagesc(normr(sortedeeg)); hold on
            plot(round(sortedrt * data.fsample),1:nTrls,'k','linewidth',1.4); hold on
            title([schar ': ' data.label{bestChan}])
            cmap = colours(100,'plasma');
            colormap(cmap);
            drawnow

            % Select just this channel and the trials of interest, in the response-locked data
            cfg = [];
            cfg.channel = RL.label{bestChan};
            cfg.trials = find(idx);
            data = ft_selectdata(cfg,RL);

            % collate response-locked data
            for trl = 1:nTrls

                y = data.trial{trl};
                x = data.time{trl};

                xidx = x >= (-T.RT(fidx(trl))-.1) & x <= 0.1;

                % baseline-correct
                if strcmp(baselineType,'resp')
                    thisbc = mean(y(:,xidx & x >= 0),2);
                    if ~any(x >= 0)
                        thisbc = y(:,xidx);
                        thisbc = mean(y(:,end));
                    end
                    y = y - thisbc;
                end

                yidx = zeros(1,length(x));
                yidx((length(yidx)-sum(xidx)+1):length(yidx)) = 1;
                yidx = yidx==1;

                % insert into X and Y variables
                if isempty(X)
                    X = nan(1,length(x));
                    X(1,yidx) = x(xidx);
                else
                    X(end+1,:) = nan(1,length(x));
                    X(end,yidx) = x(xidx);
                end

                if isempty(Y)
                    Y = nan(1,size(y,2));
                    Y(1,yidx) = y(:,xidx);
                else
                    Y(end+1,:) = nan(size(Y,2),size(Y,3));
                    Y(end,yidx) = y(:,xidx);
                end
            end
            
        % Average over central/parietal channels with positive amplitude
        case 'average'
            
            % find channels where average ERP has positive trend 1 second before RT
            thesechannels = {'C1','C2','Cz','CP1','CP2','CPz','P1','P2','Pz'};
            thisxidx = avdata.time>=-1&avdata.time<=0;
            chantrends = nan(length(thesechannels),1);
            cmap = colours(length(thesechannels),'viridis');
            figure
            for chan = 1:length(thesechannels)
                
                thischanidx = ismember(avdata.label,thesechannels{chan});
                
                winsize = 0.5; % in ms
                thiswin = [-1 (-1+winsize)];
                slidefit = [];
                while true
                    tmpxidx = avdata.time>=thiswin(1)&avdata.time<=thiswin(2);
                    thistrend = polyfit(avdata.time(tmpxidx),avdata.avg(thischanidx,tmpxidx),1);
                    slidefit = [slidefit; thistrend(1)];
                    if find(tmpxidx,1,'last') >= find(thisxidx,1,'last')
                        break
                    end
                    thiswin = avdata.time(find(avdata.time>=thiswin(1),1,'first')+1);
                    thiswin = [thiswin thiswin+winsize];
                end
                thistrend = sum(slidefit);
%                 thistrend = polyfit(avdata.time(thisxidx),avdata.avg(thischanidx,thisxidx),1);

                chantrends(chan,1) = thistrend(1);
                if thistrend(1) > 0
                    linestyle = '-';
                else
                    linestyle = ':';
                end
                
                plot(avdata.time(thisxidx),avdata.avg(thischanidx,thisxidx),'color',cmap(chan,:),'linestyle',linestyle,'linewidth',1.2); hold on
            
            end
            
            % average over positive-trending channels
            cfg = [];
            cfg.channel = thesechannels(chantrends>0);
            cfg.avgoverchan = 'yes';
            thisdata = ft_selectdata(cfg,thisdata);
            
            avdata = ft_timelockanalysis([],thisdata);
            
            plot(avdata.time(thisxidx),avdata.avg(:,thisxidx),'color','k','linewidth',2);
            title(schar)
            drawnow
            
            % organise into X and Y group variables
            channelselection.Channel{s} = cfg.channel;
            
            for trl = 1:nTrls

                y = thisdata.trial{trl};
                x = thisdata.time{trl};

                xidx = x >= (-T.RT(fidx(trl))-.1) & x <= 0.1;

                % baseline-correct
                if strcmp(baselineType,'resp')
                    thisbc = mean(y(:,xidx & x >= 0),2);
                    if ~any(x >= 0)
                        thisbc = y(:,xidx);
                        thisbc = mean(y(:,end));
                    end
                    y = y - thisbc;
                end

                yidx = zeros(1,length(x));
                yidx((length(yidx)-sum(xidx)+1):length(yidx)) = 1;
                yidx = yidx==1;

                % insert into X and Y variables
                if isempty(X)
                    X = nan(1,length(x));
                    X(1,yidx) = x(xidx);
                else
                    X(end+1,:) = nan(1,length(x));
                    X(end,yidx) = x(xidx);
                end

                if isempty(Y)
                    Y = nan(1,size(y,2));
                    Y(1,yidx) = y(:,xidx);
                else
                    Y(end+1,:) = nan(size(Y,2),size(Y,3));
                    Y(end,yidx) = y(:,xidx);
                end
            end            
    end
    
    thisxyinfo = array2table([repmat(subjects(s),nTrls,1) fidx T.Condition(idx) T.RT(idx) nan(nTrls,1)],...
        'variablenames',{'Subject','Trial','Condition','RT','StandardType'});
    thisxyinfo.StandardType = T.Type(idx);
    xyinfo = [xyinfo; thisxyinfo];

end

if strcmp(datatype,'RL') && strcmp(datamethod,'GLM')
    
    % reorganise X, Y (and add xyinfo) into concatenated data matrix
    tlengths = [];
    for s = 1:N
        tlengths(s,:) = X{s}([1 end]);
    end
    mintime = max(tlengths(:,1));
    maxtime = min(tlengths(:,2));
    
    oldx = X;
    oldy = Y;
    
    X = [];
    Y = [];
    for s = 1:N
        thisidx = oldx{s} >= mintime & oldx{s} <= maxtime;
        if isempty(X)
            X = oldx{s}(thisidx);
        else
            X = [X; oldx{s}(thisidx)];
        end
        if isempty(Y)
            Y = oldy{s}(:,thisidx);
        else
            Y = [Y; oldy{s}(:,thisidx)];
        end
    end
end

% Save
switch datamethod
    case 'GLM'
        save(fullfile(dir_save,['cppdata_' datamethod '_' datatype '_' st 'Standards_bc' upper(baselineType(1)) baselineType(2:end) '.mat']),'X','Y','xyinfo','rtbetas','betaplot','varexplained');
    case {'single','average'}
        save(fullfile(dir_save,['cppdata_' datamethod '_' datatype '_' st 'Standards_bc' upper(baselineType(1)) baselineType(2:end) '.mat']),'X','Y','xyinfo','collated','channelselection');
end

%% Plot RT Correlation

% Remove subjects with correlation < r = 0.10
[sorted,sortidx] = sort(channelselection.R);

excludeSubjects = subjects(sortidx(sorted < .1));
excludeIdx = ismember(xyinfo(:,1),excludeSubjects);

X = X(~excludeIdx,:);
Y = Y(~excludeIdx,:);
xyinfo = xyinfo(~excludeIdx,:);
collated.eeg = collated.eeg(~excludeIdx,:);
collated.rt = collated.rt(~excludeIdx,:);

% Plot
Fs = 1024;

[sorted,sortidx] = sort(collated.rt*Fs);
sortedeeg = collated.eeg(sortidx,:);
sortedeeg = sortedeeg - mean(sortedeeg(:,[1:round(0.1*Fs)]),2);
sortedeeg = normr(sortedeeg);

figure
imagesc(imgaussfilt(sortedeeg,2)); hold on
colormap(colours(100,'viridis'))
plot(sorted,1:size(sortedeeg,1),'k','linewidth',1.4); hold on
caxis(caxis/2);
set(gca,'ticklength',[0 0])

%% Plot RT

figure; histogram(xyinfo(:,4))
set(gca,'ticklength',[0 0])
xlabel('Response time (seconds)')
ylabel('Count')
title('RTs pooled across all subjects & trials')
ax = gca;
hold on
plot(repmat(median(xyinfo(:,4)),2,1),ax.YLim,'k','linewidth',1.5)
plot(repmat(quantile(xyinfo(:,4),[.25]),2,1),ax.YLim,'k--','linewidth',1.5)
plot(repmat(quantile(xyinfo(:,4),[.75]),2,1),ax.YLim,'k--','linewidth',1.5)
plot(repmat(quantile(xyinfo(:,4),[.05]),2,1),ax.YLim,'k:','linewidth',1.5)
plot(repmat(quantile(xyinfo(:,4),[.95]),2,1),ax.YLim,'k:','linewidth',1.5)
plot(repmat(quantile(xyinfo(:,4),[.01]),2,1),ax.YLim,'r:','linewidth',1.5)
plot(repmat(quantile(xyinfo(:,4),[.99]),2,1),ax.YLim,'r:','linewidth',1.5)

%% Plot beta topography

avplot = betaplot{1};
for s = 2:N
    avplot.avg(:,s) = betaplot{s}.avg;
end
avplot.avg = normr(avplot.avg')';
avplot.avg = mean(avplot.avg,2);

figure
cfg = [];
cfg.layout = 'biosemi64.lay';
ft_topoplotER(cfg,avplot);
colormap(colours(100,'viridis'))

%% SVD on pooled data

trialbytrial = true; % true to include individual trials, false to do the subject condition averages

% ------------------------------------------------------------------------------

% qrange = [.05 .95]; % interquartile range to include
% 
% q = quantile(xyinfo(:,4),qrange);
% qidx = xyinfo(:,4) >= q(1) & xyinfo(:,4) <= q(2);
% 
% allx = max(X(qidx,:));
% thisx = allx > -min(q);
% thisy = Y(qidx,thisx)';

xthreshold = 1; % need to have at least this much % of trials in each included timpoint

allx = max(X);
thisx = mean(~isnan(X)) >= xthreshold;
thisy = Y(:,thisx)';
nTrls = size(thisy,2);

zeropoint = findMin(0,allx(thisx));

% cut out trials with missing data
missingtrials = mean(isnan(thisy))>0;
disp([num2str(sum(missingtrials)) ' trials removed for time (' num2str(round((1-(nTrls/size(xyinfo,1)))*100,2)) '%)'])

thisy = thisy(:,~missingtrials);
labels = xyinfo(~missingtrials,:);
nTrls = size(labels,1);

% scale each subject's data
scaley = thisy;
for s = 1:N
    thisvec = scaley(:,labels.Subject==subjects(s));
    thismax = prctile(abs(thisvec(:)),95);
    scaley(:,labels.Subject==subjects(s)) = thisvec / thismax;
end
thisy = scaley;

% timelock SL to -1 to 0 seconds before response
newy = [];
excludeIdx = zeros(size(thisy,2),1);
for i = 1:size(thisy,2)
     thisidx = thisx & allx >= xyinfo.RT(i)-1 & allx <= xyinfo.RT(i);
     if isempty(newy)
         newy = thisy(thisidx,i);
         newy = [newy, nan(size(newy,1),size(thisy,2)-1)];
     else
         if sum(thisidx) ~= size(newy,1)
             excludeIdx(i,1) = 1;
         else
            newy(:,i) = thisy(thisidx,i);
         end
     end
end

thisy = newy(:,~excludeIdx);
labels = labels(~excludeIdx,:);
nTrls = sum(~excludeIdx);
allx = linspace(-1,0,size(thisy,1));
thisx = ones(1,length(allx))==1;

% % Average data per subject
% if ~trialbytrial
%     tmpy = [];
%     tmplabels = [];
%     cc = 0;
%     for s = 1:N
%         for c = 1:4
%             cc = cc+1;
%             thisidx = labels.Subject==subjects(s) & labels.Condition==c;
%             tmpy(cc,:) = spm_robust_average(thisy(:,thisidx),2);
%             tmplabels = [tmplabels; array2table([subjects(s) NaN c mean(labels.RT) NaN],'variablenames',labels.Properties.VariableNames)];
%         end
%     end
%     thisy = tmpy';
%     labels = tmplabels;
% end

% % Exclude trials based on robust averaging
% if trialbytrial
%     
%     % get weights from robust averaging
%     [~,W] = spm_robust_average(thisy, 2);
% 
%     % number of trials excluded, according to weighting threshold
%     weightings = 0.1:0.1:1; % what proportion of a trial (across time) has to be weighted 'bad' to be excluded entirely
%     figure
%     subject_counts = nan(N,length(weightings));
%     for w = 1:length(weightings)
%         outliers = mean(W) < weightings(w);
%         for s = 1:N
%         	subject_counts(s,w) = mean(outliers(:,labels.Subject==subjects(s)));
%         end
%         subplot(2,ceil(length(weightings)/2),w)
%         bar(subject_counts(:,w));
%         ylim([0 1])
%         title(['w = ' num2str(weightings(w))])
%     end
%     
%     % see which subjects have consistently high numbers of rejections across different thresholds
%     excludeSubjects = subjects(mean(zscore(subject_counts),2) > 2);
%     disp(['Excluding subjects [' num2str(excludeSubjects) '] based on extreme betas...'])
%     
%     % exclude these subjects
%     outliers = ismember(labels.Subject,excludeSubjects);
% 
% %     outliers = mean(W) < .8;
%     
%     thisy = thisy(:,~outliers);
%     labels = labels(~outliers,:);
%     nTrls = size(thisy,2);
% end

% Run SVD
[pred,bf,B] = svdyeahyouknowme(thisy,2);

if size(B,1)==3
    BR = B(2,:)./B(1,:);
elseif size(B,1)==4
    BR = B(3,:)./B(2,:)./B(1,:);
end

% remove trials with extreme beta ratios
q = quantile(BR,[.025 .975]);
extremetrials = BR < q(1) | BR > q(2);

if ~trialbytrial
   subj = unique(labels(extremetrials,1)); 
   extremetrials = ismember(labels(:,1),subj);
end

disp([num2str(sum(extremetrials)) ' trials removed for beta ratio (' num2str(round(mean(extremetrials)*100,2)) '%)'])

thisy = thisy(:,~extremetrials);
pred = pred(:,~extremetrials);
labels = labels(~extremetrials,:);
B = B(:,~extremetrials);
BR = BR(:,~extremetrials);

nTrls = size(thisy,2);

% assess fit of each trial
% modelfit = normr([thisy' pred']);
% modelfit = mean((modelfit(:,1:size(thisy,1)) - modelfit(:,size(thisy,1)+1:end)).^2);
modelfit = mean((thisy-pred).^2)';

% extremetrials = zscore(modelfit) > 1;
% disp([num2str(sum(extremetrials)) ' trials removed for MSE (' num2str(round(mean(extremetrials)*100,2)) '%)'])
% 
% thisy = thisy(:,~extremetrials);
% pred = pred(:,~extremetrials);
% labels = labels(~extremetrials,:);
% B = B(:,~extremetrials);
% BR = BR(:,~extremetrials);
% modelfit = modelfit(~extremetrials,:);
% nTrls = size(thisy,2);

% % Plot
% nBins = 5;
% cmap = colours(nBins,'viridis');
% figure
% for i = 1:size(B,1)-1
%     [Q,E] = discretize(B(i,:),nBins);
%     subplot(1,size(B,1)-1,i)
%     for b = 1:nBins
%         plot(allx(thisx),mean(pred(:,Q==b),2),'color',cmap(b,:),'linewidth',1.4); hold on
%     end
%     L = legend(strsplit(num2str(E(1:end-1))));
%     L.Title.String = 'Betas';
%     L.Location = 'northwest';
%     title(['Mode ' num2str(i)])
% end

cmap = [0 232 255;
        51 110 255;
        255 186 48;
        255 0 0]/255;
    
figure
for c = 1:4
    
    subplot(2,3,1)
    tmpy = mean(thisy(:,labels.Condition==c),2);
%     tmpy = tmpy - mean(tmpy);
    plot(allx(thisx),tmpy,'color',cmap(c,:),'linewidth',1.4); hold on
    title('Raw data')
    
    subplot(2,3,2)
    plot(allx(thisx),mean(pred(:,labels.Condition==c),2),'color',cmap(c,:),'linewidth',1.4); hold on
    title('Predicted data')
    
    subplot(2,3,4)
    m = mean(B(1,labels.Condition==c));
    sem = std(B(1,labels.Condition==c))/sqrt(sum(labels.Condition==c));
    plot([c c],[m+sem m-sem],'k'); hold on
    scatter(c,m,'markerfacecolor',cmap(c,:),'markeredgealpha',0); hold on
    xlim([0 5])
    title('Beta 1')
    
    subplot(2,3,5)
    m = mean(B(2,labels.Condition==c));
    sem = std(B(2,labels.Condition==c))/sqrt(sum(labels.Condition==c));
    plot([c c],[m+sem m-sem],'k'); hold on
    scatter(c,m,'markerfacecolor',cmap(c,:),'markeredgealpha',0); hold on
    xlim([0 5])
    title('Beta 2')
    
    subplot(2,3,6)
    if size(B,1) == 3
        m = mean(BR(1,labels.Condition==c));
        sem = std(BR(1,labels.Condition==c))/sqrt(sum(labels.Condition==c));
        title('Beta Ratio')
    elseif size(B,1) == 4
        m = mean(B(3,labels.Condition==c));
        sem = std(B(3,labels.Condition==c))/sqrt(sum(labels.Condition==c));
        title('Beta 3')
    end
    plot([c c],[m+sem m-sem],'k'); hold on
    scatter(c,m,'markerfacecolor',cmap(c,:),'markeredgealpha',0); hold on
    xlim([0 5])
end

% % Get slopes per trial
% slopes = nan(nTrls,1);
% intercepts = nan(nTrls,1);
% rampstart = nan(nTrls,1);
% for trl = 1:nTrls
%     
%     disp([num2str(round(trl/nTrls,2)*100) ' %'])
%     
%     smoothed = smoothdata(pred(:,trl),'loess',size(pred,1)*.3);
%     
%     tmpx = allx(thisx);
%     thistrend = polyfit(tmpx,smoothed',1);
%     if thistrend(1) > 0
%         [~,locs] = findpeaks(-smoothed);
%         if isempty(locs) && smoothed(1) < smoothed(end)
%             locs = find(smoothed==min(smoothed));
%         end
%     else
%         [~,locs] = findpeaks(smoothed);
%         if isempty(locs) && smoothed(1) > smoothed(end)
%             locs = find(smoothed==max(smoothed));
%         end
%     end
%     locs = locs(end);
%     rampstart(trl,1) = locs;
%     
%     thistrend = polyfit(tmpx(locs(end):end),smoothed(locs:end,:)',1);
%     slopes(trl,1) = thistrend(1);
%     intercepts(trl,1) = thistrend(2);
% end
% 
% figure
% for c = 1:4
%     
%     subplot(1,3,1)
%     m = mean(slopes(labels(:,3)==c));
%     sem = std(slopes(labels(:,3)==c))/sqrt(sum(labels(:,3)==c));
%     plot([c c],[m+sem m-sem],'k'); hold on
%     scatter(c,m,'markerfacecolor',cmap(c,:),'markeredgealpha',0); hold on
%     xlim([0 5])
%     title('Slopes') 
%     
%     subplot(1,3,2)
%     m = mean(intercepts(labels(:,3)==c));
%     sem = std(intercepts(labels(:,3)==c))/sqrt(sum(labels(:,3)==c));
%     plot([c c],[m+sem m-sem],'k'); hold on
%     scatter(c,m,'markerfacecolor',cmap(c,:),'markeredgealpha',0); hold on
%     xlim([0 5])
%     title('Intercepts') 
%     
%     subplot(1,3,3)
%     m = mean(rampstart(labels(:,3)==c));
%     sem = std(rampstart(labels(:,3)==c))/sqrt(sum(labels(:,3)==c));
%     plot([c c],[m+sem m-sem],'k'); hold on
%     scatter(c,m,'markerfacecolor',cmap(c,:),'markeredgealpha',0); hold on
%     xlim([0 5])
%     title('Ramp Start') 
%     
% end

% Make table
load('D:\bCFS_EEG_Reanalysis\data\Exp2\behav\stai.mat');

betas = array2table([table2array(labels(:,1:4)) modelfit B(1:2,:)' BR' thisy(zeropoint,:)' nan(nTrls,1)],...
    'variablenames',{'Subject','Trial','Condition','RT','MSE','Beta1','Beta2','BetaRatio','PeakAmp','TraitAnxiety'});

if size(B,1) == 4 % 3 modes
    betas.Beta3 = B(3,:)';
end

for S = 1:N
    betas.TraitAnxiety(betas.Subject==subjects(S)) = stai.trait(subjects(S));
end

betas.Emotion = cell(size(betas,1),1);
betas.Emotion(betas.Condition==1 | betas.Condition==2) = {'neutral'};
betas.Emotion(betas.Condition==3 | betas.Condition==4) = {'fearful'};

betas.Expectation = cell(size(betas,1),1);
betas.Expectation(betas.Condition==1 | betas.Condition==3) = {'expected'};
betas.Expectation(betas.Condition==2 | betas.Condition==4) = {'unexpected'};

betas.Emotion = categorical(betas.Emotion,unique(betas.Emotion),unique(betas.Emotion));
betas.Expectation = categorical(betas.Expectation,unique(betas.Expectation),unique(betas.Expectation));
betas.Subject = categorical(betas.Subject,unique(betas.Subject),cellstr(num2str(unique(betas.Subject))));

% % Save table
writetable(betas,fullfile(dir_save,['cppresults_' datamethod '_' datatype '_' st 'Standards_bc' upper(baselineType(1)) baselineType(2:end) '.csv']));

lme = fitlme(betas,'Beta1 ~ Emotion*Expectation + (1|Subject)')
lme = fitlme(betas,'Beta2 ~ Emotion*Expectation + (1|Subject)')
if size(B,1)==4
    lme = fitlme(betas,'Beta3 ~ Emotion*Expectation + (1|Subject)')
end
lme = fitlme(betas,'BetaRatio ~ Emotion*Expectation + (1|Subject)')

end