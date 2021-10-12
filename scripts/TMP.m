%% Get principal component for each trial, per subject (then do slope/peaktime analysis)

st = 3; % last standards
bc = 1; % baseline correct = yes
nComponents = 3; % number of basis functions for SVD

allslopes   = []; % slope leading up to response
allpeaks    = []; % 1) time of maximal negative deflection, 2) time of trough preceding response, 3) time of peak next to response
allbf       = []; % temporal basis functions (up to nComponents)
allpred     = []; % predicted data (including channel information)
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

    % Load data
    [SL,RL,T,idx] = clean_data(fullfile(dir_data,'4_epoched',[schar '_FT_epoched_' standardType{st} 'Standards_bc' num2str(baselineCorrect(bc)) '.mat'])); 

    fidx = find(idx);
    nTrls = length(fidx);
    chanidx = ~ismember(RL.label,{'M1','M2','Nz','SNz','LH','RH','LV','UV','Status'}); % contains(RL.label,{'O','P','C'}) & ~contains(RL.label,{'T','F'}); 
    nChan = sum(chanidx);
    rts = T.RT(fidx);
    conditions = T.Condition(fidx);

    % Select trials & channels
    cfg = [];
    cfg.trials = fidx;
    cfg.channel = SL.label(chanidx);
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 4;
    SL = ft_preprocessing(cfg,SL);

    % Put into data matrix
    x = SL.time{1};
    y = nan(nTrls,nChan,length(x));
    for trl = 1:nTrls
        y(trl,:,:) = SL.trial{trl};
    end
    origy = y;

    % Do SVD on all data (individual trials, NOT grouped by condition)
    thisy = reshape(y,size(y,1)*size(y,2),size(y,3))';
    [bf,~,yhat] = svdyeahyouknowme(thisy,nComponents,false);
    yhat = reshape(yhat',size(y,1),size(y,2),size(y,3));

    allbf(s,:,:) = bf';
    
    for c = 1:4
        allpred(s,c,:,:) = squeeze(mean(yhat(conditions==c,:,:)));
    end
    
    y = yhat;

    % Assess loading of each channel per trial
	betas = nan(nChan,nTrls);
    for chan = 1:nChan
        for trl = 1:nTrls
            thisy = squeeze(y(trl,chan,:));
            B = pinv([bf ones(size(bf,1),1)]) * thisy;
            betas(chan,trl) = B(1);
        end
    end
    [~,maxchan] = max(mean(betas,2));
    
    % Get peaks & slope
    peaks = nan(nTrls,4); % 1) min overall, 2) max overall, 3) min near response, 4) max near response
    slopes = nan(nTrls,2); % 1) in 500 ms before, 2) between min and max near response
    for trl = 1:nTrls
        
        tidx = findMin(0,RL.time{trl});
        
        if tidx > round(0.5*RL.fsample) % need to have at least 500 ms of pre-response data
            thisy = squeeze(y(trl,maxchan,:));

            % check sign
            y1 = normalise(bf(:,1));
            ypos = normalise(thisy);
            yneg = normalise(thisy*(-1));

            if mean(abs(y1-yneg)) < mean(abs(y1-ypos))
                thisy = thisy * (-1);
            end

            thisy = thisy(1:tidx);

            [posamp,pospks] = findpeaks(thisy);
            [negamp,negpks] = findpeaks(-thisy);

            [maxamp,maxlat] = max(thisy);
            [minamp,minlat] = max(-thisy);

            if isempty(pospks)
                [posamp,pospks] = max(thisy);
            end
            if isempty(negpks)
                [negamp,negpks] = max(-thisy);
            end

            if maxlat > minlat
                peaks(trl,1:2) = [minlat maxlat];
            end

            pkpair = [negpks(end) pospks(find(pospks > negpks(end),1,'first'))];
            if length(pkpair) > 1 && pkpair(1) < pkpair(2)
                peaks(trl,3:4) = pkpair; 
            end

            if ~all(isnan(peaks(trl,3:4)))
                l = polyfit(pkpair(1):pkpair(2),thisy(pkpair(1):pkpair(2))',1);
                slopes(trl,1) = l(1);
            else
                slopes(trl,1) = NaN;
            end

            thisx = (tidx - round(0.5*RL.fsample)):tidx; % 500 ms prior to response
            l = polyfit(thisx,thisy(thisx)',1);
            slopes(trl,2) = l(1);
        end
    end
    
    for c = 1:4
        allpeaks(s,c,:) = nanmean(peaks(conditions==c,:));
        allslopes(s,c,:) = nanmean(peaks(conditions==c,:));
    end
end

% plot basis set
figure
for i = 2:nComponents
    subplot(1,nComponents-1,i-1)
    y = squeeze(mean(allbf(ismember(subjects,theseSubjects),:,:)))';
    y = y(:,1) + y(:,i)*linspace(-2^i,2^i,10);
%     y = y + abs(y(1,:));
    for j = 1:size(y,2)
        y(:,j) = normalise(y(:,j));
    end
    P = plot(x,y,'linewidth',1.4);
    cmap = colours(length(P),'viridis');
    for i = 1:length(P)
        P(i).Color = cmap(i,:);
    end
    set(gca,'ticklength',[0 0])
    xlim(x([1 end]))
end

figure
cmap = [0 232 255;
        51 110 255;
        255 186 48;
        255 0 0]/255;
        
subplot(1,2,1)
for c = 1:4
    y = squeeze(allpred(ismember(subjects,theseSubjects),c,:));
%     y = y-mean(y(:,x<0),2); % set the start to 0
    m = mean(y);
    sem = std(y)/sqrt(thisN);
    upper = m+sem;
    lower = m-sem;
%     patch([x fliplr(x)],[upper fliplr(lower)],cmap(c,:),'facealpha',.2,'edgecolor','none'); hold on
    plot(x,m,'color',cmap(c,:),'linewidth',1.5); hold on
    set(gca,'ticklength',[0 0])
end
xlim(x([1 end]))

rtx = squeeze(nanmean(rtlog));
rtx(2,:) = rtx(2,:)/max(rtx(2,:));
tidx = rtx(1,:) > -1;
subplot(1,2,2)
for c = 1:4
    y = squeeze(allrl(ismember(subjects,theseSubjects),c,tidx));
%     y = y-mean(y(:,x<0),2); % set the start to 0
    m = mean(y);
    sem = std(y)/sqrt(thisN);
    upper = m+sem;
    lower = m-sem;
%     patch([x fliplr(x)],[upper fliplr(lower)],cmap(c,:),'facealpha',.2,'edgecolor','none'); hold on
    plot(rtx(1,tidx),m,'color',cmap(c,:),'linewidth',1.5); hold on
    set(gca,'ticklength',[0 0])
end

% Peak times
figure
for i = 1:3
    subplot(1,3,i)
    y = zscore(squeeze(allpeaks(:,:,i))')';
    for c = 1:4
        m = mean(y(ismember(subjects,theseSubjects),c));
        sem = std(y(ismember(subjects,theseSubjects),c))/sqrt(thisN);
        upper = m+sem;
        lower = m-sem;
        scatter(c,m); hold on
        plot([c c],[upper lower],'k'); hold on
    end
    xlim([0 5])
    title(['Peak ' num2str(i)])
end

% Slopes
figure
for i = 1:2
    subplot(1,2,i)
    y = zscore(squeeze(allslopes(:,:,i))')';
    for c = 1:4
        m = mean(y(ismember(subjects,theseSubjects),c));
        sem = std(y(ismember(subjects,theseSubjects),c))/sqrt(thisN);
        upper = m+sem;
        lower = m-sem;
        scatter(c,m); hold on
        plot([c c],[upper lower],'k'); hold on
    end
    xlim([0 5])
    title(['Peak ' num2str(i)])
end