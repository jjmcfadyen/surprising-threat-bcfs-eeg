d = readtable('D:\bCFS_EEG_Reanalysis\results\exp3.csv');
anx = d.Anxiety;
zanx = zscore(anx);
d = table2array(d(:,1:4));

anxmean = 39.49701;

%% Anxiety scatter plots to overlay LME trendlines

low = zanx <= -0.5;
medium = zanx > -0.5 & zanx < 0.5;
high = zanx >= 0.5;

figure
for i = 1:3

    if i==1
        idx = low;
        thisrange = linspace(min(anx)+anxmean, anx(findMin(-0.5,zanx))+anxmean,3);
        xc = thisrange(2);
    elseif i ==2
        idx = medium;
        xc = anx(findMin(0,zanx))+anxmean;
    elseif i==3
        idx = high;
        thisrange = linspace(anx(findMin(0.5,zanx))+anxmean,max(anx)+anxmean,3);
        xc = thisrange(2);
    end

    x = linspace(xc-3,xc+3,4);
%     draw_boxpointviolin(d(idx,:),'xtick',x,'width',1,'cmap',cmap)
    for c = 1:4
        m = mean(d(idx,c));
        sem = std(d(idx,c))/sqrt(sum(idx));
        plot(repmat(x(c),2,1),[m+sem m-sem],'k','linewidth',1.2); hold on
        scatter(x(c),m,60,'markerfacecolor',cmap(c,:),'markeredgecolor','none'); hold on
    end

end

set(gca,'ticklength',[0 0])
xlim([20 60])

%% Anxiety correlations

figure
for i = 1:2

    if i==1
        thiscmap = cmap(1:2,:);
    elseif i==2
        thiscmap = cmap(3:4,:);
    end

    x = anx + anxmean;
    y = table2array(d(:,5+i));
    lm = polyfit(x,y,1);
    lm = x*lm(1) + lm(2);

    subplot(1,2,i)
    scatter(x(y>0),y(y>0),40,'markerfacecolor',thiscmap(1,:),'markerfacealpha',.5,'markeredgecolor','k','linewidth',1); hold on
    scatter(x(y<0),y(y<0),40,'markerfacecolor',thiscmap(2,:),'markerfacealpha',.5,'markeredgecolor','k','linewidth',1); hold on
    plot(x,lm,'k','linewidth',2); hold on
    set(gca,'ticklength',[0 0])
    xlim([10 70])
    ylim([-0.4 0.8])
    plot([20 60],[0 0],'k:','linewidth',1.5); hold on

end
