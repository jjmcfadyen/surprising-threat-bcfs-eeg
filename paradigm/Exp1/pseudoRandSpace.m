%% function for pseudo-randomising spacing of target trials across block (or stimuli within rsvp stream)

% numTargets:   number of targets across block
% numTrials:    number of trials in block
% minGap:       minimum number of distractor trials between targets


function [targetOrder, adjSpaces] = pseudoRandSpace(numTargets,numTrials,minGap,plotGaps)

% variables
if nargin<3 || isempty(minGap)
    minGap = 0;
end
if nargin<4
    plotGaps = 0;
end

avgSpaceIncludingTarget = numTrials/numTargets;
avgSpaceBetweenTargets = avgSpaceIncludingTarget-1;
maxGap = avgSpaceBetweenTargets-minGap;

% pseudo-randomise how much to adjust spaces between targets
adjSpaces = abs(normrnd(0, maxGap/2, 1, numTargets)); % one-tailed random gaussian distribution
adjSpaces = max(0, adjSpaces); % remove adjustments that are too small
adjSpaces = min(maxGap, adjSpaces); % remove adjustments that are too large

% allocate targets
targetOrder = false(1,numTrials);
for targ = 1:numTargets
    centredLoc = 1 + (targ-1) * avgSpaceIncludingTarget;
    adjustedLoc = centredLoc + adjSpaces(targ);
    targetOrder(round(adjustedLoc)) = 1;
end

% check code has done what it should, if not do again
if length(targetOrder)~=numTrials || ...              % number of trials matches input
    sum(targetOrder)~=numTargets   ;                    % number of targets matches input
    disp(targetOrder)
    disp(numTargets)
    disp(numTrials)
    warning([num2str(length(targetOrder)-numTrials) ' less targets than requested!!!'])
end
if any(diff(find(targetOrder))<=minGap) % check minimum gap has been applied
    error('check code')
end


% plot gap distribution?
if plotGaps
    gaps = diff(find(targetOrder))-1;
    for t = 1:max(gaps)+1
        test(t) = sum(gaps==(t-1));
    end
    figure
    bar(0:max(gaps),test)
    xlabel('gaps between targets')
end






%% SUPERCEDED



% % linear distribution
% linearDist = linspace(adjRange(1),adjRange(2),numTargets);
% adjSpaces = linearDist(randperm(numTargets)); 



% % work out number of trials between each target
% if strcmp(distribution,'gaussian') % gaussian distribution of target locations
% 
%     % determine target spacing using gaussian distribution
%     while ~ok
%         spacing = normrnd(avgSpace,avgSpace/2,1,numTargets); % gaussian distribution of spaces between targets
%         spacing = max(0,spacing); % remove <=0
%         spacing = min(avgSpace*2,spacing); % remove >=avgSpace*2
%         spacing = round(spacing)+1; % round and add 1 for target location itself
%         ok = sum(spacing)<=numTrials; % check if number of spaces (and targets) will fit into number of trials
%     end
%     
% elseif isnumeric(spacingRange) && isequal(size(spacingRange),[1 2]) % linear distribution of target locations
%     
%     % determine target spacingu using linear distribution
%     while ~ok
%         spacing = pseudoRand(range(spacingRange)+1,numTargets) + spacingRange(1)-1; % rand(1,numTargets)*diff(spacingRange)+spacingRange(1);
%         spacing = spacing + 1; % add 1 for target location itself
%         ok = sum(spacing)<=numTrials; % check if number of spaces (and targets) will fit into number of trials
%     end
%     
% else 
%     error('check range format')
% end


