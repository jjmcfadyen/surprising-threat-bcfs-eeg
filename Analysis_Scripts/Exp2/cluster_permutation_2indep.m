function stat = cluster_permutation_2indep(data1,data2,neighbours,timewindow,channel,these_subjects)

    % Cluster-based permutation test
    cfg = [];
    cfg.neighbours = neighbours;
    cfg.latency = timewindow;
    cfg.channel = channel; % 'EEG'
    cfg.parameter = 'avg';
    cfg.method = 'montecarlo';
    cfg.statistic = 'indepsamplesT';
    
    Nsub = [length(these_subjects) length(these_subjects)];
    cfg.design(1,:) = [ones(1,Nsub(1)), ones(1,Nsub(2))*2];
    cfg.design(2,:) = [1:Nsub(1) 1:Nsub(2)];
    cfg.ivar = 1;

    cfg.correctm = 'no';
    cfg.clusteralpha = .05;
    cfg.correcttail = 'prob';
    cfg.numrandomization = 1000;
    
    stat = ft_timelockstatistics(cfg,data1{these_subjects},data2{these_subjects});
    
end
