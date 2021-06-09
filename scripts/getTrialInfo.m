function [trialinfo,nTrls] = getTrialInfo(D)

nTrls = length(D.events);
events = D.events;    

vNames = {'trialStart','trialEnd','responseTime','stimulusType','responseType'};
trialinfo = array2table(nan(nTrls,length(vNames)),'variablenames',vNames);

parfor trl = 1:nTrls

    thisEvent = events{trl};
    thisEvent = thisEvent(~contains(extractfield(thisEvent,'type'),'artefact_OSL')); % ignore artefact event markers

    onsetEvent    = thisEvent(ismember(extractfield(thisEvent,'value'),1:4)); % start of trial
    offsetEvent   = thisEvent(extractfield(thisEvent,'value') >= 99); % end of trial
    if isempty(offsetEvent)
        offsetEvent = [];
        offsetEvent.time = onsetEvent.time+3;
    end

    responseEvent = thisEvent(ismember(extractfield(thisEvent,'value'),21:24)); % response
    if isempty(responseEvent)
        responseEvent = [];
        responseEvent.time = NaN;
        responseEvent.value = NaN;
    end

    % Fill table
    thistable = array2table(nan(1,length(vNames)),'variablenames',vNames);
    
    thistable.trialStart(1)     = onsetEvent.time;
    thistable.trialEnd(1)       = offsetEvent.time;
    thistable.responseTime(1)   = responseEvent.time;
    thistable.stimulusType(1)   = onsetEvent.value;
    thistable.responseType(1)   = responseEvent.value-20;
    
    trialinfo(trl,:) = thistable;
    
end

% Mark responses
trialinfo.missed = isnan(trialinfo.responseType);
trialinfo.RT = trialinfo.responseTime - trialinfo.trialStart;

% Mark bad trials
trialinfo.badTrials = zeros(nTrls,1);
trialinfo.badTrials(D.badtrials) = 1;

trialinfo.badTrials(trialinfo.RT < .5) = 1; % mark responses faster than 500 ms as artefacts

% Label standards as different types (first, middle, or last)
trialinfo.standardType = cell(nTrls,1);
for trl = 1:nTrls
    if ismember(trialinfo.stimulusType(trl),[2 4])
        trialinfo.standardType{trl} = 'deviant';
    else
        if trl==1 % first trial is first standard
            trialinfo.standardType{trl} = 'first';
        elseif trl==size(trialinfo,1) % last trial is last standard
            trialinfo.standardType{trl} = 'last';
        else
            if ismember(trialinfo.stimulusType(trl-1),[2 4]) % if trial before was a deviant
                trialinfo.standardType{trl} = 'first';
            elseif ismember(trialinfo.stimulusType(trl+1),[2 4]) % if trial after is a deviant
                trialinfo.standardType{trl} = 'last';
            else
                trialinfo.standardType{trl} = 'middle';
            end
        end
    end
end

end