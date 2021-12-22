function plot_ddm()

ddm = readtable('D:\bCFS_EEG_Reanalysis\results\ddmparams.csv'); 

cmap = [0 232 255;
        51 110 255;
        255 186 48;
        255 0 0]/255;

figure
for i = 1:3
    
    if i==1
        thisy = ddm(:,contains(ddm.Properties.VariableNames,{'nondecision'}));
    elseif i==2
        thisy = ddm(:,contains(ddm.Properties.VariableNames,{'drift'}));
    elseif i==3
        thisy = ddm(:,contains(ddm.Properties.VariableNames,{'boundary'}));
    end
    thisy = table2array(thisy);

    subplot(1,3,i)
    for c = 1:4

        binwidth = diff(linspace(min(thisy(c,:)),max(thisy(c,:)),20));
        binwidth = binwidth(1);

        [x,y] = beeswarm(thisy(:,c),binwidth,0.15);

        q = quantile(y,[.25 .75]);
        patch([c-.15 c-.15 c+.15 c+.15],[q(1) q(2) q(2) q(1)],cmap(c,:),'facealpha',.25,'edgecolor',cmap(c,:)); hold on
        plot([c c],[q(2) max(thisy(:,c))],'color',cmap(c,:),'linewidth',1.2); hold on
        plot([c c],[q(1) min(thisy(:,c))],'color',cmap(c,:),'linewidth',1.2); hold on

%         ds = normalise(ksdensity(y))*0.35;
%         patch([ds -fliplr(ds)]+c,[linspace(min(y),max(y),length(ds)) fliplr(linspace(min(y),max(y),length(ds)))],cmap(c,:),'facealpha',.25,'edgecolor',cmap(c,:),'edgealpha',.7); hold on
%         scatter(x+c,y,20,'markerfacecolor',cmap(c,:),'markeredgecolor','none','markerfacealpha',.5); hold on

        m = mean(y);
        sem = std(y)/sqrt(size(y,1));
        upper = m+sem;
        lower = m-sem;

        plot([c c],[lower upper],'k','linewidth',1.3); hold on
        scatter(c,m,70,'markerfacecolor',cmap(c,:),'markeredgecolor','k'); hold on

    end
    xlim([0 5])
    ylim([-2 2])
    set(gca,'ticklength',[0 0])

     if i==1
        title('Nondecision')
    elseif i==2
        title('Drift')
    elseif i==3
        title('Boundary')
    end

end

end