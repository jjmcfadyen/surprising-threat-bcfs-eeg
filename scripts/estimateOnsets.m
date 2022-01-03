function [trlonsets rltrlonsets] = estimateOnsets(Y,X,rt,conditions,nIterations)
% estimateOnsets(Y,X)
% Estimates two types of onsets: the first (i.e., the mode of onset times up until response) and last (i.e., the most significant onset up until response)
% Variables:
%       Y               = data (trials x channels x samples)
%       X               = time, in seconds (1 x samples)
%       rt              = response times (trials x 1)
%       conditions      = conditions (trials x 1)
%       nIterations     = number of iterations to randomly permute response time (set to 0 to ignore)

%% Set parameters

onsettimemin = 0.1; % can't have an onset before 100 ms
onsetlatewindow = 0.5; % how many seconds before RT to look for motor/decision-related increase in activity

%% Get info

nTrls = size(Y,1);
nChan = size(Y,2);

%% Estimate onsets

% Crop data
xidx = X>onsettimemin & X<3;
Y = Y(:,:,xidx);
X = X(xidx);

% Estimate onsets per trial
trlonsets = nan(nTrls,nChan,2,nIterations+1);
rltrlonsets = nan(nTrls,nChan,2,nIterations+1);
parfor trl = 1:nTrls
    disp(['Trial ' num2str(trl)])
    thisrt = rt(trl);
    for chan = 1:nChan

%         % estimate onsets for an ever-increasing window (e.g., 0 to 1, 0 to 10, 0 to 20, 0 to 30, ... 0 to 3)
%         thisonset = [];
%         tp = X(2);
%         tplist = [tp];
%         while round(tp,3) <= thisrt
% 
%             thisdata = squeeze(Y(trl,chan,X<=tp));
%             o = findchangepts(thisdata,'statistic','linear');
%             if isempty(o)
%                 o = NaN;
%             end
%             thisonset = [thisonset o];
% 
%             try
%                 tp = X(findMin(tp,X)+10);
%                 if tp <= thisrt
%                     tplist = [tplist tp];
%                 end
%             catch
%                 tp = Inf;
%             end
%         end
% 
%         tmp_trlonsets = nan(2,nIterations+1);
%         tmp_rltrlonsets = nan(2,nIterations+1);
%         for i = 1:nIterations+1
% 
%             if i==1
%                 thisrt = rt(trl);
%             else
%                 thisrt = randsample(0.9:0.01:3,1); %rt(randsample(setdiff(1:length(rt),trl),1)); % either take a random RT from 0.9 to 3 seconds, or randomly pick another trial's RT
%             end
% 
%             lasto = thisonset(findMin(thisrt,tplist)); % find most significant build up in activity in a window from 0 to response
% %             firsto = thisonset(findMin(tplist,thisx(lasto)));
% %             firsto = mode(thisonset(tplist <= thisx(lasto)));
% 
%             if ~isnan(lasto)
%                 lasto = thisx(lasto);
%             end
%             if ~isnan(firsto)
%                 firsto = thisx(firsto);
%             end
%     
%             tmp_trlonsets(:,i) = [firsto lasto]';
%             tmp_rltrlonsets(:,i) = -(thisrt - [firsto lasto])';
%         end
% 
%         % estimate "last" onset as most significant change in 1 second window preceding response
%         xidx = X>=(thisrt-onsetlatewindow) & X<thisrt;
%         thisx = X(xidx);
%         lasto = thisx(findchangepts(squeeze(Y(trl,chan,xidx)),'statistic','linear'));
% 
%         % estimate "first" onset as most significant change prior to the last most significant change
%         xidx = X>onsettimemin & X<lasto;
%         thisx = X(xidx);
%         firsto = thisx(findchangepts(squeeze(Y(trl,chan,xidx)),'statistic','linear'));

        xidx = X>onsettimemin & X<(thisrt);
        thisx = X(xidx);
        o = findchangepts(squeeze(Y(trl,chan,xidx)),'statistic','linear','maxnumchanges',3);
        if length(o) < 3
            o = findchangepts(squeeze(Y(trl,chan,xidx)),'statistic','linear','maxnumchanges',2);
        end
        if length(o) >= 2
            firsto = thisx(o(1));
            lasto = thisx(o(end));
        else
            firsto = NaN;
            lasto = NaN;
        end

        tmp_trlonsets = [firsto lasto];
        tmp_rltrlonsets = -(thisrt - [firsto lasto])';

        trlonsets(trl,chan,:,:) = tmp_trlonsets;
        rltrlonsets(trl,chan,:,:) = tmp_rltrlonsets;
    end
end
for i = 1:2
    if i==1
        disp('::: FIRST =====================================')
    elseif i==2
        disp('::: LAST =====================================')
    end
    tmp1 = squeeze(trlonsets(:,:,i,1));
    tmp2 = squeeze(rltrlonsets(:,:,i,1));
    disp([':::::::::::: ' num2str(round(mean(mean(isnan(tmp1)))*100,2)) '% of trials with no identifiable onset time'])
    disp([':::::::::::: Mean onset time = ' num2str(round(nanmean(tmp1(:)),3)) ' seconds'])
    disp([':::::::::::: Mean onset time = ' num2str(round(nanmean(tmp2(:)),3)) ' seconds response locked)'])
    disp([':::::::::::: Mean no. of onsets after response = ' num2str(nanmean(tmp2(:)>=0))])
end

end