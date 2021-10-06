function [SL,RL,T,idx] = clean_data(filename)
% Takes the preprocessed data, interpolates bad channels, and sorts into stimulus-locked and response-locked time series

addpath('D:\Toolboxes\fieldtrip-20191119')
ft_defaults

cfg = [];
cfg.method = 'template';
cfg.layout = 'biosemi64.lay';
neighbours = ft_prepare_neighbours(cfg);

% Load data
tmp = load(filename);
D = tmp.thisD;
T = tmp.thisT;
idx = tmp.idx;
clear tmp

if size(T,1) ~= length(D.trial)
    error('Mismatch with behavioural file')
end

% Check for bad channels & interpolate
[badchannels,SL] = fix_badchannels(D,neighbours);

% Shift time axis for response-locked data
offset = nan(length(SL.trial),1);
for trl = 1:length(SL.trial)
    offset(trl,1) = -findMin(T.RT(trl),SL.time{trl});
end

cfg = [];
cfg.offset = offset;
RL = ft_redefinetrial(cfg,SL);

end