function matched = matchBehaviour(behav,meg)
% behav = trial information from 'preprocess_behav.m' for ONE subject
% meg   = trial info from 'getTrialInfo.m' for ONE subject

B = [behav.Condition, round(behav.RT,3)]; % condition number (1 to 4), response time
M = [meg.stimulusType, round(meg.RT,3)]; % condition number (1 to 4), response time

nM = size(meg,1);
meg.behavIdx = nan(nM,1);

% Match on condition code & RT
conMatch = M(:,1) == B(:,1);
rtMatch = M(:,2) > B(:,2)-.1 & M(:,2) < B(:,2)+.1;

matched = true;
if any(~conMatch)
    matched = false;
end

end