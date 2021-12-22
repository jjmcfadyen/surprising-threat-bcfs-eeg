function draw_boxpointviolin(data,varargin)
% draw_boxpointviolin(data,varargin)
% Draws a boxplot on the left, and a dot plot & violin on the right of an x point
%
% data      = data matrix (rows = observations, columns = conditions/blocks/whatever)
% width     = how wide the box & violin plot extends (default = 0.2)
% xtick     = points along the x axis for each column in 'data' to fall (default = 1:size(data,2))
% cmap      = color map for each column in x (if only one color, will be applied to all)
% drawlines = draw lines between scatter plots along x axis (default = false)
% linecolor = colour of lines (if drawlines is true) (default = black, or cmap if cmap has only one colour defined)

% Requires 'beeswarm.m' and 'normalise.m' and 'colours.m' functions (usually found in 'utils' directory)

%% Get parameters

nCol = size(data,2);

% defaults
width = 0.2;
xtick = 1:nCol;
if nCol<3
    cmap = repmat([0 0 0],nCol,1); 
else
    cmap = colours(nCol,'viridis');
end
markersize = 20;
drawlines = false;
linecolor = [];

for v = 1:length(varargin)
    if ischar(varargin{v})
        switch varargin{v}
            case 'width'
                width = varargin{v+1};
            case 'xtick'
                xtick = varargin{v+1};
            case 'cmap'
                cmap = varargin{v+1};
            case 'markersize'
                markersize = varargin{v+1};
            case 'linecolor'
                linecolor = varargin{v+1};
            case 'drawlines'
                drawlines = varargin{v+1};
        end
    end
end

if size(cmap,1)==1 & nCol>1
    cmap = repmat(cmap,nCol,1);
end

if size(unique(cmap,'rows'),1)==1
    linecolor = cmap(1,:);
elseif isempty(linecolor)
    linecolor = [0 0 0];
end

%% Draw

pairs = [];
for b = 1:nCol

    % scatter points & 1/2 violin & boxplot
    binwidth = diff(linspace(min(data(:)),max(data(:)),20));
    binwidth = binwidth(1);

    [x,y] = beeswarm(data(:,b),binwidth,width/2);

    ds = normalise(ksdensity(y))*width;
    patch([ds ds(1)]+xtick(b),[linspace(min(y),max(y),length(ds)) min(y)],cmap(b,:),'facealpha',.25,'edgecolor',cmap(b,:),'edgealpha',.7); hold on

    scatter(x+xtick(b)+(width/2),y,markersize,'markerfacecolor',cmap(b,:),'markeredgecolor',cmap(b,:),'markerfacealpha',.5,'markeredgealpha',.2); hold on
    pairs = [pairs (x+xtick(b)+(width/2)) y];

    q = quantile(y,[.25 .75]);
    patch([xtick(b)-width xtick(b)-width xtick(b) xtick(b) xtick(b)-width],...
        [q(1) q(2) q(2) q(1) q(1)],cmap(b,:),'facealpha',.3,'edgecolor',cmap(b,:),'edgealpha',.7); hold on
    plot([xtick(b)-width xtick(b)],repmat(nanmean(y),2,1),'color',cmap(b,:),'linewidth',1.5); hold on
    plot(repmat(xtick(b)-(width/2),2,1),quantile(y,[.05,.95]),'color',cmap(b,:),'linewidth',1.8); hold on
end

if drawlines
    for k = 1:size(pairs,1)
        plot(pairs(k,[1:2:end]),pairs(k,[2:2:end]),'color',[linecolor .2],'linewidth',1.2); hold on
    end
end
end