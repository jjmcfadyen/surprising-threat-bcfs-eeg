function [badchannels,interpolated] = fix_badchannels(D,neighbours)

% get data info
channels = find(~ismember(D.label,{'M1','M2','Nz','SNz','LH','RH','LV','UV','Status'}));
    
nTrls = length(D.trial);
nChan = length(channels);

% organise data into 3D matrix
d = nan(nChan,length(D.time{1}),nTrls);
for trl = 1:nTrls
    d(:,:,trl) = D.trial{trl}(channels,:);
end

% use GESD to determine channel outliers
d = std(reshape(d,size(d,1),size(d,2)*size(d,3)),[],2); % get standard deviation of each channel

outliers = gesd(d,.005); % relatively more strict alpha criterion for what makes a channel an outlier

badchannels = D.label(channels(outliers));

% interpolate
cfg = [];
cfg.method = 'weighted';
cfg.badchannel = badchannels;
cfg.neighbours = neighbours;
interpolated = ft_channelrepair(cfg,D);

end