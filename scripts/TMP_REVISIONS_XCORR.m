st=3;
b=2;

theseSubjects = setdiff(subjects,[3]); % subject 3 missed about half the trials
thisN = length(theseSubjects);

grand = [];
grandz = [];
for s = 1:thisN

    subject = subjects(theseSubjects(s));
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

    % Organise data
    X = nan(length(SL.trial),length(SL.label),length(SL.time{1}));
    for trl = 1:length(SL.trial)
        X(trl,:,:) = SL.trial{trl};
    end

    X = X(idx,:,:);
    Y = T.Condition(idx);

    chanidx = ~ismember(SL.label,{'Nz','SNz','LH','RH','LV','UV','Status','M1','M2'});
    X = X(:,chanidx,:);

    % Save averages to variable
    m = [];
    for c = 1:4
        m(c,:,:) = mean(X(Y==c,:,:));
    end

    grand(s,:,:,:) = m;

    % Scale each channel/trial and save to variable
    xidx = SL.time{1} >= 0;
    Z = [];
    for trl = 1:size(X,1)
        for chan = 1:size(X,2)
            Z(trl,chan,:) = zscore(X(trl,chan,xidx));
        end
    end

    m = [];
    for c = 1:4
        m(c,:,:) = mean(Z(Y==c,:,:));
    end

    grandz(s,:,:,:) = m;

end

%% Cross correlate

% average across channels
channels = SL.label(chanidx);

% y = squeeze(mean(grand(:,:,chanidx,:),3));
y = squeeze(grand(:,:,ismember(channels,'CPz'),:));
m = squeeze(mean(y));
sem = squeeze(std(y)/sqrt(size(y,1)));

cmap = [77, 242, 255;
        98, 174, 255;
        255, 180, 29;
        255, 68, 142]/255;

x = SL.time{1};
figure
for c = 1:4
    patch([x fliplr(x)],[m(c,:)+sem(c,:) fliplr(m(c,:)-sem(c,:))],cmap(c,:),'facealpha',.2,'edgecolor','none'); hold on
    plot(x,m(c,:),'color',cmap(c,:),'linewidth',1.5); hold on
end

figure
for s = 1:size(y,1)
    subplot(5,6,s)
    for c = 1:4
        plot(x,squeeze(y(s,c,:)),'color',cmap(c,:),'linewidth',1.5); hold on
    end
end





% cross correlate
y = squeeze(mean(grand(:,:,ismember(SL.label,{'Cz','CPz'}),:),3));

xc = [];
for s = 1:size(y,1)

    y1 = squeeze(y(s,1,:));
    y2 = squeeze(y(s,2,:));
    [C,lags] = xcorr(y1,y2);
    xc(s,1) = lags(C==max(C));

    y1 = squeeze(y(s,3,:));
    y2 = squeeze(y(s,4,:));
    [C,lags] = xcorr(y1,y2);
    xc(s,2) = lags(C==max(C));
end

figure
for chan = 1:size(grand,3)
    subplot(8,8,chan)
    draw_boxpointviolin(squeeze(xc(:,chan,:)));
    [h,p] = ttest(squeeze(xc(:,chan,1)),squeeze(xc(:,chan,2)));
    title([SL.label{chan} ' p=' num2str(round(p,3))])
end


pvals = [];
for chan = 1:size(crosscorr,2)
    [h,p] = ttest(squeeze(crosscorr(:,chan,1)),squeeze(crosscorr(:,chan,2)));
    pvals(chan,1) = p;
end

%% Strategies

output = preprocess_behav(1,false);
T = output.trialdata;
T.Experiment = repmat(1,size(T,1),1);

output = preprocess_behav(2,false);
output.trialdata.Subject = output.trialdata.Subject + max(T.Subject);
output.trialdata.Experiment = repmat(2,size(output.trialdata,1),1);
T = [T; output.trialdata];

% check accuracy
tsubjects = unique(T.Subject);
tN = length(tsubjects);
acc = nan(tN,4);
missed = nan(tN,4);
for s = 1:tN
    for c = 1:4
        acc(s,c) = mean(T.Acc(T.Subject==tsubjects(s) & T.Condition==c));
        missed(s,c) = mean(isnan(T.RT(T.Subject==tsubjects(s) & T.Condition==c)));
    end
end

expvec = [];
for s = 1:tN
    expvec(s,1) = unique(T.Experiment(T.Subject==tsubjects(s)));
end

acc = array2table([acc expvec],'variablenames',{'EN','UN','EF','UF','Experiment'});
writetable(acc,'acc.csv')

missed = array2table([missed  expvec],'variablenames',{'EN','UN','EF','UF','Experiment'});
writetable(acc,'missed.csv')




accpvals = [];
for c1 = 1:4
    for c2 = 1:4
        [~,accpvals(c1,c2)] = ttest(acc(:,c1),acc(:,c2));
    end
end

missedpvals = [];
for c1 = 1:4
    for c2 = 1:4
        [~,missedpvals(c1,c2)] = ttest(missed(:,c1),missed(:,c2));
    end
end

figure
subplot(1,2,1)
draw_boxpointviolin(acc);
subplot(1,2,2)
draw_boxpointviolin(missed);
