function [LI,LO,nFolds] = bc_getFolds(Y,foldsize)
% Mini function to evenly sample trials into different folds

conditions = unique(Y);
nCon = length(conditions);
nTrials = size(Y,1);


LI = {}; % left in
LO = {}; % left out

% first fold is all data included
LI{1} = 1:nTrials;

% evenly sample from conditions until all trials have been used
selectedtrials = zeros(nTrials,1); % log how many times each trial is included in a fold
f = 1;
while true
    
    f = f + 1;
    
    LO{f} = nan(nCon,foldsize);
    for c = 1:nCon
        thisrange = find(Y==conditions(c) & selectedtrials==min(selectedtrials(Y==conditions(c))));
        if length(thisrange) < foldsize
            thesetrials = thisrange;
            thisrange = setdiff(find(Y==conditions(c) & selectedtrials<=min(selectedtrials(Y==conditions(c)))+1),thisrange);
            thesetrials = [thesetrials; datasample(thisrange,foldsize-length(thesetrials),'replace',false)];
        else
            thesetrials = datasample(thisrange,foldsize,'replace',false);
        end
        LO{f}(c,:) = thesetrials;
        selectedtrials(thesetrials) = selectedtrials(thesetrials) + 1;
    end
    LI{f} = setdiff(1:nTrials,LO{f}(:));
    
    if ~any(selectedtrials==0)
        break
    end
end
nFolds = f;

end