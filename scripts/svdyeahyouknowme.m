function [BF,VF,Yhat] = svdyeahyouknowme(Y,nBF,makePlots)
% nBF = Number of basis functions
% Y = time x space

%==========================================================================
% Derive basis set and apply to data
%==========================================================================
oM = mean(Y);
mcY = Y - oM;

[u,s,v]     = svd(mcY);                 % SVD
ds          = diag(s);                  % Get scales for basis set
BF          = u(:,1:nBF).*ds(1:nBF)';   % Apply scales to components of choice i.e up to nBF (temporal basis function)
VF          = v(1:nBF,:)'*ds(1:nBF);    % spatial basis function

% if betas trend negatively, then flip them around
bftrend = polyfit(1:size(BF,1),BF(:,1)',1);
if bftrend(1) < 0
    BF = BF * (-1);
end

bf          = [BF ones(length(BF),1)];    % Add constant term for GLM
B           = pinv(bf) * Y;               % GLM
Yhat        = bf(:,1:end-1)*B(1:end-1,:); % Reconstructed responses
Yhat        = Yhat + oM;                  % put the mean back in

%==========================================================================
%                       *******I Love you********
%==========================================================================

if makePlots
    figure
    
    subplot(2,2,1)
    plot(bf(:,1:end-1),'linewidth',1.4)
    title('Basis Set')
    xlabel('Time')
    ylabel('a.u.')
    axis square
    
    subplot(2,2,2)
    plot(ds/sum(ds),'k','linewidth',2)
    ylabel('Variance Explained')
    xlim([0 7])
    xticks(1:7)
    axis square
    title('Eigenvalue')
    
    subplot(2,2,3)
    hold on
    P = plot(bf(:,1) + bf(:,2)*linspace(-1,1,10),'linewidth',1.4);
    cmap = colours(length(P),'viridis');
    for i = 1:length(P)
        P(i).Color = cmap(i,:);
    end
    plot(bf(:,1),'k','linewidth',2.2)
    title('First basis scaled by second basis')
    axis square
    
    if nBF > 2
        subplot(2,2,4)
        hold on
        P = plot(bf(:,1) + bf(:,3)*linspace(-1,1,10),'linewidth',1.4);
        cmap = colours(length(P),'viridis');
        for i = 1:length(P)
            P(i).Color = cmap(i,:);
        end
        plot(bf(:,1),'k','linewidth',1.2)
        axis square
        title('First basis scaled by third basis')
    else
        subplot(2,2,4)
        plot(mean(Y,2))
        title('Raw data (average)')
    end
end

% figure
% random_sample = randsample(1:size(Y,2),36,false);
% for i = 1:36
%     subplot(6,6,i)
%     hold on
%     plot(Y(:,random_sample(i)),'k','linewidth',1.4)
%     plot(Yhat(:,random_sample(i)),'r','linewidth',1.4)
%     title(num2str(B(1,random_sample(i)) - B(end,random_sample(i))))
%     axis square
% end
% suptitle('Random sample of 36 trials reconstructed')
% legend({'Raw','GLM'})
%
%
% figure
%
% subplot(1,2,1)
%
% x = 1:size(Y,1);
% m = mean(Y,2)';
% sem = std(Y')/2;%/sqrt(length(x));
% upper = m+sem;
% lower = m-sem;
%
% patch([x fliplr(x)],[upper fliplr(lower)]','k','facealpha',.2,'edgealpha',0); hold on
% plot(m,'k','linewidth',1.6); hold on
% title('Raw')
% axis square
%
%
% subplot(1,2,2)
%
% x = 1:size(Yhat,1);
% m = mean(Yhat,2)';
% sem = std(Yhat');%/sqrt(length(x));
% upper = m+sem;
% lower = m-sem;
%
% patch([x fliplr(x)],[upper fliplr(lower)]','r','facealpha',.2,'edgealpha',0); hold on
% plot(m,'r','linewidth',1.6); hold on
% axis square
% title('GLM reconstructed')
% suptitle('All trials')

end
