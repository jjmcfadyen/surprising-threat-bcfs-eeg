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
onsetrtbuffer = 0.3; % don't look at activity in the 100 ms leading up to response
stepsize = 0.05; % step size in increasing window, in ms

%% Get info

nTrls = size(Y,1);
nChan = size(Y,2);

%% Estimate onsets

X = round(X,3);

% Crop data
xidx = X>=0 & X<3;
Y = Y(:,:,xidx);
X = X(xidx);

% Estimate onsets per trial
trlonsets = nan(nTrls,nChan,2,nIterations+1);
rltrlonsets = nan(nTrls,nChan,2,nIterations+1);
parfor trl = 1:nTrls

    disp(['Trial ' num2str(trl)])
    thisrt = rt(trl);
    thisidx = X>=onsettimemin & X<=(thisrt-onsetrtbuffer);
    thisx = X(thisidx);
    thisy = squeeze(Y(trl,:,thisidx));

    for chan = 1:nChan

        % find most significant change using ever-increasing window
        fwdtp = thisx(1)+stepsize;
        bwdtp = thisx(end)-stepsize;
        fwdonsets = [];
        bwdonsets = [];
        fwdtplist = [];
        bwdtplist = [];
        while true
            
            fwdtplist = [fwdtplist fwdtp];
            bwdtplist = [bwdtplist bwdtp];

            o = findchangepts(thisy(chan,thisx<=fwdtp),'statistic','linear');
            if isempty(o)
                o = NaN;
            else
                tmpx = thisx(thisx<=fwdtp);
                o = tmpx(o);
            end
            fwdonsets = [fwdonsets o];

            o = findchangepts(thisy(chan,thisx>=bwdtp),'statistic','linear');
            if isempty(o)
                o = NaN;
            else
                tmpx = thisx(thisx>=bwdtp);
                o = tmpx(o);
            end
            bwdonsets = [bwdonsets o];

            fwdtp = fwdtp + stepsize;
            bwdtp = bwdtp - stepsize;

            if fwdtp > max(thisx) || bwdtp < min(thisx)
                break
            end

        end

%         % looking backwards (biaseded towards response time), find the point at which the change times can be split into two groups
%         thisidx = ~isnan(bwdonsets);
%         tmpx = bwdtplist(thisidx);
%         lasto = tmpx(findchangepts(bwdonsets(thisidx),'statistic','mean'));
% 
%         % looking forwards (biased towards stimulus onset), find the point at which the change times can be split into two groups BEFORE the response-biased change time
%         thisidx = ~isnan(fwdonsets) & fwdtplist < lasto;
%         tmpx = fwdtplist(thisidx);
%         firsto = tmpx(findchangepts(fwdonsets(thisidx),'statistic','mean'));

        % take the most frequent onst time in each vector
        [N, ~, BIN] = histcounts(fwdonsets,'binwidth',stepsize);
        [~,maxbin] = max(N);
        firsto = mean(fwdonsets(BIN==maxbin));

        [N, ~, BIN] = histcounts(bwdonsets,'binwidth',stepsize);
        [~,maxbin] = max(N);
        lasto = mean(bwdonsets(BIN==maxbin));

        if firsto > lasto
            
            tmp = fwdonsets(fwdtplist<firsto);
            [N, ~, BIN] = histcounts(tmp,'binwidth',stepsize);
            [~,maxbin] = max(N);
            firsto = mean(tmp(BIN==maxbin));

            tmp = bwdonsets(bwdtplist>lasto);
            [N, ~, BIN] = histcounts(tmp,'binwidth',stepsize);
            [~,maxbin] = max(N);
            lasto = mean(tmp(BIN==maxbin));

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
    disp([':::::::::::: Mean onset time = ' num2str(round(nanmean(tmp2(:)),3)) ' seconds response locked'])
    disp([':::::::::::: Mean no. of onsets after response = ' num2str(nanmean(tmp2(:)>=0))])
end
end