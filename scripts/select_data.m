function [SL,RL,T,idx] = select_data(filename,neighbours,standardType,dataType)
% [SL,RL,T,idx] = select_data(filename,neighbours,standardType,dataType)
% Takes the preprocessed data, interpolates bad channels, and sorts into stimulus-locked and response-locked time series

% addpath('D:\Toolboxes\fieldtrip-20191119')
% ft_defaults

%%

% Load data
tmp = load(filename);
switch dataType
    case 'orig'
        SL = tmp.thisD;
    case 'scd_finite'
        SL = tmp.scd_finite;
    case 'scd_spline'
        SL = tmp.scd_spline;
    case 'scd_hjorth'
        SL = tmp.scd_hjorth;
end
T = tmp.thisT;
trialinfo = tmp.trialinfo;
nTrls = length(SL.trial);
clear tmp

if size(T,1) ~= length(SL.trial)
    error('Mismatch with behavioural file')
end

% Put data into 3D matrix
d = nan(length(SL.trial),length(SL.label),length(SL.time{1}));
for trl = 1:nTrls
    d(trl,:,:) = SL.trial{trl};
    d(trl,:,SL.time{trl} > T.RT(trl)) = NaN; % ignore data after response
end

% Z-score all data
mu = nanmean(d(:));
sigma = nanstd(d(:));

d = nan(length(SL.trial),length(SL.label),length(SL.time{1}));
for trl = 1:nTrls
    SL.trial{trl} = (SL.trial{trl} - mu) / sigma;
    d(trl,:,:) = SL.trial{trl};
    d(trl,:,SL.time{trl} > T.RT(trl)) = NaN; % ignore data after response
end

% Shift time axis for response-locked data
offset = nan(nTrls,1);
for trl = 1:nTrls
    offset(trl,1) = -findMin(T.RT(trl),SL.time{trl});
end

cfg = [];
cfg.offset = offset;
RL = ft_redefinetrial(cfg,SL);

%% Get index of trials matching this standard type (all, first, last)

% Only select trials not marked as bad
idx = trialinfo.badTrial==0;

% Only select trials that have a response
idx = idx & ~isnan(T.RT); % right response triggers missing, so use behavioural log instead

% Only select trials with correct responses
idx = idx & T.Acc;

% Remove trials where RTs are too fast (< 500 ms)
idx = idx & T.RT > .5;

% Remove trials where the RT was > 3 SDs from the mean RT
rtoutliers = find(idx);
rtoutliers(:,2) = abs(zscore(T.RT(idx))) > 3;
idx(rtoutliers(rtoutliers(:,2)==1)) = 0;

% Select the type of standard
switch standardType
    case 'first'
        idx = idx & (contains(trialinfo.standardType,'deviant') | contains(trialinfo.standardType,'first'));
    case 'last'
        idx = idx & (contains(trialinfo.standardType,'deviant') | contains(trialinfo.standardType,'last'));
end

% use GESD to determine trial outliers
sd = nanstd(reshape(d,size(d,1),size(d,2)*size(d,3)),[],2); % get standard deviation of each trial
outliers = gesd(sd,.05); % relatively strict alpha criterion for what makes a trial an outlier
disp(['Removing ' num2str(sum(outliers)) ' outliers...'])

idx = idx & ~outliers';

end