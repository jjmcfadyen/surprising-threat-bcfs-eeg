% Visualise the observed vs. simulated data according to specific model

function vis_sim(rt,sim_rt)

    % Convert to hazard function
    Q = quantile(rt,[.001 .97]);
    x = Q(1):.001:Q(2);

    try
        fit = fitdist(rt(rt >= Q(1) & rt <= Q(2)),'Burr');
    catch
        fit = fitdist(rt,'Burr');
    end
    hazard = pdf('Burr',x,fit.alpha,fit.c,fit.k)./(1-cdf('Burr',x,fit.alpha,fit.c,fit.k));

    try
        sim_fit = fitdist(sim_rt(sim_rt >= Q(1) & sim_rt <= Q(2)),'Burr');
        sim_hazard = pdf('Burr',x,sim_fit.alpha,sim_fit.c,sim_fit.k)./(1-cdf('Burr',x,sim_fit.alpha,sim_fit.c,sim_fit.k));
    catch sim_hazard = nan(length(x),1);
    end

    % Convert to quantiles
    Q = [.05:.05:.95];
    Q_rt = quantile(rt,Q);
    Q_sim = quantile(sim_rt,Q);

    % plot
    figure;
    subplot(1,3,1); histogram(rt); hold on; histogram(sim_rt); xlabel('RT'); ylabel('Frequency'); legend({'Observed','Simulated'})
    subplot(1,3,2); plot(x,hazard); hold on; plot(x,sim_hazard); ylabel('Failure Time (s)'); xlabel('Hazard Rate'); legend({'Observed','Simulated'})
    subplot(1,3,3); scatter(Q,Q_rt); hold on; scatter(Q,Q_sim);
    plot(Q,Q_rt,'k'); plot(Q,Q_sim,'k');
    xlabel('Response Proportion'); ylabel('RT Quantiles (s)'); legend({'Observed','Simulated'})

    disp(['Mean RT = ' num2str(nanmean(rt))])
    disp(['Mean Simulated RT = ' num2str(nanmean(sim_rt))])

end
