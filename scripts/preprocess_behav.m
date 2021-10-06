function output = preprocess_behav(expnum,makeplots)

%% Set parameters

% Directories
dir_home = 'D:\bCFS_EEG_Reanalysis\';
addpath(fullfile(dir_home,'scripts'));
addpath(genpath(fullfile(dir_home,'scripts')));

dir_data = fullfile(dir_home,'data',['Exp' num2str(expnum)],'behav');

cmap = [77, 242, 255;
        98, 174, 255;
        255, 180, 29;
        255, 68, 142]/255;
    
conditionLabels = {'EN','UN','EF','UF'}';
conditionCodes  = [31 32
                   21 22
                   11 12
                   41 42];
nCon = length(conditionLabels);

load(fullfile(dir_data,'stai.mat')); % loads 'stai' variable
               
%% Load data

load(fullfile(dir_data,'stai.mat'));

filelist = dir(fullfile(dir_data,'s*trial_info.mat'));

subjects = cell(length(filelist),1);
for f = 1:length(filelist)
   tmp = strsplit(filelist(f).name,'_');
   schar = tmp{1}(2:end);
   if length(schar) == 1
       schar = ['S0' schar(end)];
   else
       schar = ['S' schar];
   end
   subjects{f,1} = schar;
end
N = length(subjects);

trialdata = [];
conditiondata = [];
trialcount = [];
demographics = array2table(nan(N,3),'variablenames',{'Subject','Age','Sex'});
for s = 1:N
    
    schar = [];
    if s < 10
        schar = ['S0' num2str(s)];
    else
        schar = ['S' num2str(s)];
    end
    
    % read in data
    f = find(contains(subjects,schar));
    d = load(fullfile(filelist(f).folder,filelist(f).name));
    d = d.trial_info;
    
    demographics.Subject(s) = d.subject.num;
    demographics.Age(s) = d.subject.age;
    demographics.Sex(s) = d.subject.gender;
    
    subject = d.subject.num;
    
    nTrls = numel(d.trials);
    nBlocks = size(d.trials,1);
    trialsPerBlock = size(d.trials,2);
    
    % format long data
    vN = {'Subject','Block','Trial','Condition','Emotion','Expectation','Type','RT','Orientation','Response','Acc'};
    thislong = array2table(nan(nTrls,length(vN)),'variablenames',vN);
    
    thislong.Subject = repmat(subject,nTrls,1);
    thislong.Block = sort(repmat(1:nBlocks,1,trialsPerBlock))';
    thislong.Trial = repmat(1:trialsPerBlock,1,nBlocks)';
    thislong.Condition = reshape(d.trials',nTrls,1);
    thislong.RT = reshape(d.rt',nTrls,1);
    thislong.Acc = reshape(d.acc'==1,nTrls,1);
    
    for c = 1:nCon
        thislong.Condition(ismember(thislong.Condition,[conditionCodes(c,:)])) = c;
    end
    
    thislong.Emotion = cell(nTrls,1);
    thislong.Emotion(ismember(thislong.Condition,[1:2])) = {'neutral'};
    thislong.Emotion(ismember(thislong.Condition,[3:4])) = {'fearful'};
    
    thislong.Expectation = cell(nTrls,1);
    thislong.Expectation(ismember(thislong.Condition,[1 3])) = {'expected'};
    thislong.Expectation(ismember(thislong.Condition,[2 4])) = {'unexpected'};
    
    if isfield(d,'orientation')
        thislong.Orientation = reshape(d.orientation',nTrls,1);
        thislong.Response(thislong.Acc) = thislong.Orientation(thislong.Acc);
        thislong.Response(~thislong.Acc & thislong.Orientation==1) = 2;
        thislong.Response(~thislong.Acc & thislong.Orientation==2) = 1;
    end
    
    thislong.DominantEye = repmat(d.dominanteye,nTrls,1);
    
    thislong.Outlier = zeros(nTrls,1);
    thislong.Outlier(abs(zscore(thislong.RT))>=5) = 1;
    thislong.Outlier = thislong.Outlier==1;
    
    try
        thislong.StateAnxiety   = repmat(stai.state(s),nTrls,1);
        thislong.TraitAnxiety   = repmat(stai.trait(s),nTrls,1);
        thislong.STAI           = repmat(stai.full(s),nTrls,1);
    catch
        thislong.StateAnxiety   = nan(nTrls,1);
        thislong.TraitAnxiety   = nan(nTrls,1);
        thislong.STAI           = nan(nTrls,1);
    end
    
    % label standard types (first, middle, last)
    thislong.Type = cell(nTrls,1);
    
    deviantIdx = find(thislong.Condition == 2 | thislong.Condition == 4);
    
    lastIdx = deviantIdx - 1;
    lastIdx(lastIdx <= 0,:) = [];
    
    firstIdx = deviantIdx + 1;
    firstIdx(firstIdx > nTrls,:) = [];
    
    middleIdx = setdiff(1:nTrls,[firstIdx; lastIdx; deviantIdx]);
    
    thislong.Type = cell(nTrls,1);
    thislong.Type(deviantIdx) = {'deviant'};
    thislong.Type(firstIdx) = {'first'};
    thislong.Type(lastIdx) = {'last'};
    thislong.Type(middleIdx) = {'middle'};
    
    if isempty(trialdata)
        trialdata = thislong;
    else
        trialdata = [trialdata; thislong];
    end
    
    % identify outliers by getting typical response time for correct trials
    thisidx = ~isnan(thislong.RT) & thislong.Acc & thislong.RT > .5;
    thislong.Outliers = zeros(size(thislong,1),1);
    thislong.Outliers(thisidx) = abs(zscore(thislong.RT(thisidx))) > 3;
    
    % get trial counts and data for plotting
    for i = 1:3 % standard types (all, first, last)
        for c = 1:4

            idx = thislong.Condition==c & thislong.Acc==1 & ~isnan(thislong.RT);

            if i==2
                idx = idx & (contains(thislong.Type,'deviant') | contains(thislong.Type,'first'));
            elseif i==3
                idx = idx & (contains(thislong.Type,'deviant') | contains(thislong.Type,'last'));
            end

            trialcount(subject,i,c) = sum(idx);

            for m = 1:2 % mean, median
                if m==1
                    y = mean(thislong.RT(idx));
                elseif m==2
                    y = median(thislong.RT(idx));
                end

                conditiondata(subject,m,i,c) = y;
            end
        end
    end
end

%% Plot

conditionorder = [1 4 3 2]; % EN UF, EF UN
meanormedian = 'median';

if makeplots
    
    figure
    cc = 0;
    allylims = [];
   
    switch meanormedian
        case 'mean'
            m = 1;
        case 'median'
            m = 2;
    end
    
    for i = 1:3

        cc = cc + 1;
        subplot(1,3,cc)

        for c = 1:4
            thismean = squeeze(mean(conditiondata(:,m,i,conditionorder(c))));
            thiserror = squeeze(std(conditiondata(:,m,i,conditionorder(c))))/sqrt(size(conditiondata,1));
            upper = thismean+thiserror;
            lower = thismean-thiserror;
            plot([c c],[lower upper],'k'); hold on
            scatter(c,thismean,'markerfacecolor',cmap(conditionorder(c),:),'markeredgealpha',0); hold on
        end
        xlim([0 5])

        if m==1
            ylabel('mean');
        elseif m==2
            ylabel('median');
        end

        if i==1
            title('all standards')
        elseif i==2
            title('first standard')
        elseif i==3
            title('last standards')
        end

        ax = gca;
        allylims = [allylims; ax.YLim];

        set(gca,'XTick',1:nCon);
        set(gca,'XTickLabels',conditionLabels(conditionorder))
        set(gca,'ticklength',[0 0])

    end
    
    for i = 1:cc
        subplot(1,3,cc)
        ylim([min(allylims(:,1)) max(allylims(:,2))])
    end
    
    sgtitle(['Experiment ' num2str(expnum)])
    
end

%% Store in output

output = [];
output.trialdata = trialdata;
output.conditiondata = conditiondata;
output.conditiondatadimensions = {'subject','mean or median','standard type: all, first, middle, last','condition: EN, UN, EF, UF'};
output.trialcount = trialcount;

output.demographics = demographics;

if expnum==1 % last subject is missing the STAI info
    output.demographics.STAI = [stai.full; NaN];
    output.demographics.State = [stai.state; NaN];
    output.demographics.Trait = [stai.trait; NaN];
else
    output.demographics.STAI = stai.full;
    output.demographics.State = stai.state;
    output.demographics.Trait = stai.trait;
end

end