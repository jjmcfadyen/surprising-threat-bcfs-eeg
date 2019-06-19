clear all
close all
clc

%% Variables

experiment = 2;

% maximum trial length
if experiment == 1
    max_rt = 10;
elseif experiment == 2
    max_rt = 3;
end

parameter_names = {'a','Ter','eta','st','v'};
s = .15; % noise in random walk

%% Load data

T = readtable(['D:\Scratch\bCFS_EEG_Reanalysis\results\trial_data_exp' num2str(experiment) '.csv'],'FileType','spreadsheet','ReadVariableNames',true);

subjects = unique(T.Subject);

conditions = zeros(size(T,1),1);
conditions(T.Emotion == 1 & T.Expectation == 1) = 1;
conditions(T.Emotion == 1 & T.Expectation == 2) = 2;
conditions(T.Emotion == 2 & T.Expectation == 1) = 3;
conditions(T.Emotion == 2 & T.Expectation == 2) = 4;
nCon = length(unique(conditions));
condition_names = {'EN','UN','EF','UF'};

colours = [ 87 150 231
    120 62 223
    196 33 132
    241 134 33
    ]/256;

display_subjects = false; % set to true to see output of individual subjects

%% Group

% rt = T.RT;
% cidx = conditions;
% [sorted_rt, sort_idx] = sort(rt); % sort by RT (fastest to slowest)
% 
% % Display Quantiles
% Q = [.05:.05:.95];
% 
% figure
% trial_freq = nan(nCon,length(Q));
% for con = 1:nCon
%     
%     this_rt = rt(cidx==con);
%     Q_rt = quantile(this_rt,Q);
%     
%     % Plot
%     plot(Q,Q_rt,'Color',colours(con,:),'LineWidth',2); hold on
%     scatter(Q,Q_rt,'MarkerFaceColor',colours(con,:),'MarkerEdgeColor','k'); hold on
%     
%     % Store frequency of trials in each quantile category
%     for q = 1:length(Q)
%         trial_freq(con,q) = sum(this_rt <= Q(q));
%     end
%     
% end
% xlabel('Response Proportion')
% ylabel('RT Quantiles (s)')
% 
% % Hazard Function
% Q = quantile(rt,[.001 .97]);
% x = Q(1):.001:Q(2);
% 
% figure
% for con = 1:nCon
%     
%     this_rt = rt(cidx==con);
%     this_rt = this_rt(this_rt >= Q(1) & this_rt <= Q(2));
%     fit = fitdist(this_rt,'Burr');
%     
%     hazard = pdf('Burr',x,fit.alpha,fit.c,fit.k)./(1-cdf('Burr',x,fit.alpha,fit.c,fit.k));
%     
%     if con == 1 || con == 3
%         plot(x,hazard,'Color',colours(con,:),'LineWidth',2); hold on
%     elseif con == 2 || con == 4
%         plot(x,hazard,'Color',colours(con,:),'LineWidth',2,'LineStyle','--'); hold on
%     end
%     ylabel('Failure Time (s)')
%     xlabel('Hazard Rate')
%     
% end
% legend(condition_names);
% 
% % ESTIMATE PARAMTERS
% % Fixed parameters
% tau = .001; % time step size, in seconds
% nSamples = 20000;
% 
% a = 1;        % boundary (fixed)
% 
% % start_v = (a*(.1/std(rt))^2)^(1/3);
% % start_Ter = mean(rt) / a*start_v;
% % start_st = start_Ter*.25;
% % start_eta = start_v*.5;
% 
% group_parameter_estimates = nan(nCon,5);
% group_chi_square = nan(1,nCon);
% group_sim_data = cell(1,nCon);
% group_obs_data = cell(1,nCon);
% for con = 1:nCon
% 
%     disp('----------------------------------------------------------')
%     disp(['--- CONDITION ' num2str(con) '-------------------------------'])
%     disp('----------------------------------------------------------')
%     
%     this_rt = rt(cidx == con);
% 
% %     start_v = (a*(.1/std(rt))^2)^(1/3);
% %     start_Ter = mean(rt) / a*start_v;
% %     start_st = start_Ter*.25;
% %     start_eta = start_v*.5;
% %     start_params = [a start_Ter start_eta start_st start_v];
%     start_params = [a 1 .5 .5 1];
%     
%     minvals = [.5 0.01 0.01 0.01 0.01];
%     maxvals = [1.5 2 1 1 2];
%     
%     loss = @(x)loglikelihood(this_rt,max_rt,[x s]'); % input the parameters here
%     options = optimset('Algorithm','interior-point','TolCon',1e-15,'TolX',1e-15,'TolFun',1e-15,'FunValCheck','on','FinDiffType','central');
%     
%     [mle_params, logl] = fmincon(loss,start_params,[],[],[0 1 0 .5 0],min(this_rt)-.001,minvals,maxvals,[],options);
%     
%     mle_params(5) = mle_params(5) * a/mle_params(1);
%     mle_params(3) = mle_params(3) * a/mle_params(1);
%     mle_params(1) = a;
% 
%     % Visualise initial guess
%     sim_rt = simulateProcess([mle_params s],nSamples,max_rt,tau);
%     vis_sim_newversion(this_rt,sim_rt);
%     title(['Initial Guess for Group: ' condition_names{con}])
%     
%     % Parameter estimation
%     minvals = [a 0.01 0.01 0.01 0.01];
%     maxvals = [a 1 .5 .5 1];
%     loss = @(x)modelFit(this_rt,[x s],nSamples,tau,max_rt);
%     nTries = 4;
%     targetChiSquare = 30;
% 
%     cc = 0;
%     carryon = true;
%     log = zeros(nTries,length(minvals)+1);
%     while carryon
%         cc = cc  +1;
%         if cc >= nTries
%             carryon = false;
%         end
%         disp(['(((((((( TRY #' num2str(cc) ')))))))))'])
%         [opt_params, chiSquare] = fminsearchbnd(loss,mle_params,minvals,maxvals,...
%             optimset('TolFun',0.5,'TolX',0.05,'MaxFunEvals',150,'Display','iter'));
%         log(cc,:) = [chiSquare, opt_params];
%         if chiSquare < targetChiSquare
%             carryon = false;
%         end
%     end
%     opt_params = log(log(:,1) == min(log(:,1)),2:end);
%     opt_params(3) = opt_params(3)/a;
%     opt_params(5) = opt_params(5)/a;
%     group_chi_square(1,con) = chiSquare;
% 
%     % Visualise optimal model
%     sim_rt = simulateProcess([opt_params s],nSamples,max_rt,tau);
%     vis_sim_newversion(this_rt,sim_rt(~isnan(sim_rt)));
%     title(['Optimised Model for Group: ' condition_names{con}])
% 
%     % Save
%     group_sim_data{1,con} = sim_rt;
%     group_obs_data{1,con} = this_rt;
%     group_parameter_estimates(con,:) = opt_params;
% 
% end
% 
% % disp('--- GROUP AVERAGE: ')
% % for con = 1:nCon
% %     for p = 1:size(group_parameter_estimates,2)
% %         disp(['----------------- ' condition_names{con} ', ' parameter_names{p} ' = ' num2str(group_parameter_estimates(con,p))])
% %     end
% % end
% 
% % group_parameter_estimates(:,end+1) = group_parameter_estimates(:,5)./group_parameter_estimates(:,3);
% 
% figure
% jitter = linspace(-.2,.2,nCon);
% for p = 1:size(group_parameter_estimates,2)
%     plot(p+jitter,group_parameter_estimates(:,p),'Color',[.3 .3 .3]); hold on 
% end
% for con = 1:nCon
%     scatter((1:length(parameter_names))+jitter(con),group_parameter_estimates(con,:),'filled','MarkerFaceColor',colours(con,:)); hold on
% end
% xlim([0 size(group_parameter_estimates,2)+1])
% set(gca,'XTick',1:size(group_parameter_estimates,2))
% set(gca,'XTickLabels',parameter_names)
% xlabel('Parameter Type')
% ylabel('Estimated Value')
% title('Parameter Estimates (estimated across entire group per condition)')
% 
% save('group_estimation_v3.mat','group_parameter_estimates','group_chi_square','group_sim_data','group_obs_data');
% 
% %% Try to recover parameters
% 
% sim_parameters      = [1 2 .2 .5 1
%                        1 2 .5 .5 3
%                        1 .5 .2 .1 1
%                        1 .5 .4 .1 3];
% 
% trial_range = [60 300 700 20000];
% 
% recovery_iterations = 20;
% retrieved_params = {};
% for S = 1%:size(sim_parameters,1)
%     
%     these_sim_params = sim_parameters(S,:);
%     
%     for T = 1:length(trial_range)
%         
%         this_trial_range = trial_range(T);
%         retrieved_params{S,T} = [];
%         
%         for it = 1:recovery_iterations
% 
%             % Simulate data
%             this_rt = simulateProcess([these_sim_params s],this_trial_range,max_rt,tau);
% 
%             % Use MLE for estimating initial guess
%             start_v = (a*(.1/std(this_rt))^2)^(1/3);
%             start_Ter = mean(this_rt) / a*start_v;
%             start_st = start_Ter*.25;
%             start_eta = start_v*.5;
% 
%             start_params = [a start_Ter start_eta start_st start_v];
%             minvals = [.5 0.1 0.01 0.01 0.01];
%             maxvals = [1.5 max_rt max_rt max_rt max_rt];
% 
%             ub = start_params > maxvals;
%             lb = start_params < minvals;
%             start_params(ub) = maxvals(ub);
%             start_params(lb) = minvals(lb);
% 
%             loss = @(x)loglikelihood(this_rt,max_rt,[x s]');
%             options = optimset('Algorithm','interior-point','TolCon',1e-15,'TolX',1e-15,'TolFun',1e-15,'FunValCheck','on','FinDiffType','central');
%             [mle_params, logl] = fmincon(loss,start_params,[],[],[0 1 0 .5 0],min(this_rt)-.001,minvals,maxvals,[],options);
%             mle_params(5) = mle_params(5) * a/mle_params(1);
%             mle_params(3) = mle_params(3) * a/mle_params(1);
%             mle_params(1) = a;
% 
% %             % Visualise initial guess
% %             sim_rt = simulateProcess([mle_params s],nIterations,max_rt,tau);
% %             vis_sim_newversion(this_rt,sim_rt);
% %             title('initial guess')
% 
%             % Parameter estimation
%             minvals = [a 0.3 0.01 0.01 0.01];
%             maxvals = [a max_rt/2 max_rt/4 max_rt/4 max_rt/2];
%             loss = @(x)modelFit(this_rt,[x s],nSamples,tau,max_rt);
%             nTries = 4;
%             targetChiSquare = 30;
% 
%             cc = 0;
%             carryon = true;
%             while carryon
%                 [opt_params, chiSquare] = fminsearchbnd(loss,mle_params,minvals,maxvals,...
%                     optimset('TolFun',0.5,'TolX',0.05,'MaxFunEvals',150,'Display','iter'));
%                 cc = cc  +1;
%                 if cc > nTries || chiSquare < targetChiSquare
%                     carryon = false;
%                 end
%                 disp(['(((((((( TRY #' num2str(cc) ')))))))))'])
%             end
%             opt_params(3) = opt_params(3)/a; % this doesn't really do anything anyway (a is 1)...
%             opt_params(5) = opt_params(5)/a;
% 
% %             % Visualise optimal model
% %             sim_rt = simulateProcess([opt_params s],nIterations,max_rt,tau);
% %             vis_sim_newversion(this_rt,sim_rt);
% %             title(['Group: ' condition_names{con}])
% 
%             % Compare
%             retrieved_params{S,T} = [retrieved_params{S,T}; opt_params];
% 
%         end
%     end
% end
% 
% save('recovery_test.mat','retrieved_params','sim_parameters','trial_range');
% 
% % Plot parameter recovery
% figure
% nParams = size(starting_parameters,2);
% cc = 0;
% for S = 1:size(retrieved_params,1)
%     for T = 1:size(retrieved_params,2)
%         cc = cc + 1;
%         this_it = size(retrieved_params{S,T},1);
%         jitter = linspace(-.3,.3,this_it);
%         subplot(size(retrieved_params,1),size(retrieved_params,2),cc)
%         for p = 1:nParams
%             scatter(ones(1,this_it)*p+jitter,retrieved_params{S,T}(:,p),'filled','MarkerFaceAlpha',.8); hold on
%             scatter(p,sim_parameters(1,p),70,'filled','k','Marker','s'); hold on
%         end
%         xlim([.5 size(retrieved_params{S,T},2)+.5])
%         title(['Starting parameter set ' num2str(S) ', ' num2str(trial_range(T)) ' trials'])
%         set(gca,'XTick',1:nParams)
%         set(gca,'XTickLabels',parameter_names)
%     end
% end

%% Per subject

equate_trials = true;
nIterations = 10;

parameter_estimates = nan(length(subjects),nCon,length(parameter_names));
all_chi_square = nan(length(subjects),nCon,10);
all_sim_data = cell(length(subjects),nCon);
all_obs_data = cell(length(subjects),nCon);
for subj = 1:length(subjects)
    
    subject = subjects(subj);
    
    cidx = conditions(T.Subject == subject);
    trlidx = T.ExpTrial(T.Subject == subject);
    rt = T.RT(T.Subject == subject);
    
%     start_v = a/mean(rt);
%     start_Ter = min(rt)/2;
%     start_st = start_Ter/4;
%     start_eta = start_v/2;
% 
%     % Use MLE for estimating initial guess
%     start_params = [a start_Ter start_eta start_st start_v];
%     minvals = [.5 0.01 0.01 0.01 0.01];
%     maxvals = [1.5 max_rt/min(rt) max_rt/min(rt)/2 mean(rt)/2 mean(rt)];
    
%     if display_subjects 
%         
%         %% Display Quantiles
% 
%         Q = [.05:.05:.95];
% 
%         figure
%         trial_freq = nan(nCon,length(Q));
%         for con = 1:nCon
% 
%             this_rt = rt(cidx==con);
%             Q_rt = quantile(this_rt,Q);
% 
%             % Plot
%             plot(Q,Q_rt,'Color',colours(con,:),'LineWidth',2); hold on
%             scatter(Q,Q_rt,'MarkerFaceColor',colours(con,:),'MarkerEdgeColor','k'); hold on
% 
%             % Store frequency of trials in each quantile category
%             for q = 1:length(Q)
%                 trial_freq(con,q) = sum(this_rt <= Q(q));
%             end
% 
%         end
%         xlabel('Response Proportion')
%         ylabel('RT Quantiles (s)')
% 
%         %% Hazard Function
% 
%         Q = quantile(rt,[.001 .97]);
%         x = Q(1):.001:Q(2);
% 
%         figure
%         for con = 1:nCon
% 
%             this_rt = rt(cidx==con);
%             this_rt = this_rt(this_rt >= Q(1) & this_rt <= Q(2));
%             fit = fitdist(this_rt,'Burr');
% 
%             hazard = pdf('Burr',x,fit.alpha,fit.c,fit.k)./(1-cdf('Burr',x,fit.alpha,fit.c,fit.k));
% 
%             if con == 1 || con == 3
%                 plot(x,hazard,'Color',colours(con,:),'LineWidth',2); hold on
%             elseif con == 2 || con == 4
%                 plot(x,hazard,'Color',colours(con,:),'LineWidth',2,'LineStyle','--'); hold on
%             end
%             ylabel('Failure Time (s)')
%             xlabel('Hazard Rate')
% 
%         end
%         legend(condition_names);
%     end
                    
    %% Estimate parameters
    
    % Following: https://doi.org/10.1145/2554850.2554872 & https://doi.org/10.1111/jsr.12166
    % https://github.com/amiyapatanaik/one-choice-DDM
    
    % Fixed parameters
    tau = .001; % time step size, in seconds
    nSamples = 20000;
    a = 1;        % boundary (fixed)

    min_trials = [];
    for c = 1:nCon
        min_trials(c) = sum(cidx == c);
    end
    con_min = find(min_trials == min(min_trials));
    min_trials = min_trials(con_min(1));
    
    tic
    for con = 1:nCon
        
        con_rt = rt(cidx == con);
        
        if equate_trials == true && length(con_rt) ~= min_trials
            these_iterations = nIterations;
        else these_iterations = 1;
        end
        
        for n = 1:these_iterations
                
            initial_fit = true;
            ii = 0;
            while initial_fit
                
                ii = ii + 1;
                if ii > 100
                    warning(['Cannot get good estimate of initial fit for subject ' num2str(subject)])
                    mle_params = start_params;
                end
                
                if length(con_rt) > min_trials && equate_trials
                    keep_sampling = true;
                    cc = 0;
                    while keep_sampling
                        this_rt = con_rt(datasample(1:length(con_rt),min_trials,'replace',false));
                        try 
                            Q = quantile(this_rt,[.001 .97]);
                            fit = fitdist(this_rt(this_rt >= Q(1) & this_rt <= Q(2)),'Burr');
                            keep_sampling = false;
                        end
                        disp('Randomly sampling data for this condition until Burr distribution is satisfied...')
                        cc = cc + 1;
                        if cc > 100
                            keep_sampling = false; 
                        end
                    end
                else this_rt = con_rt;
                end

                start_v = a/mean(this_rt);
                start_Ter = min(this_rt)/2;
                start_st = start_Ter/4;
                start_eta = start_v/4;

                % Use MLE for estimating initial guess
                start_params = [a start_Ter start_eta start_st start_v];
                minvals = [.5 0.01 0.01 0.01 0.01];
%                 maxvals = [1.5 mean(this_rt)/2 min(this_rt)/2 min(this_rt)/2 mean(this_rt)];
                maxvals = [1.5 min(this_rt) a/mean(this_rt) min(this_rt)/2 1/(min(this_rt)/2)];
                ub = start_params > maxvals;
                lb = start_params < minvals;
                start_params(ub) = maxvals(ub);
                start_params(lb) = minvals(lb);

                loss = @(x)loglikelihood(this_rt,max_rt,[x s]'); % input the parameters here
                options = optimset('Algorithm','interior-point','TolCon',1e-15,'TolX',1e-15,'TolFun',1e-15,'FunValCheck','on','FinDiffType','central');

                [mle_params, logl] = fmincon(loss,start_params,[],[],[0 1 0 .5 0],min(this_rt)-.001,minvals,maxvals,[],options);
                mle_params(5) = mle_params(5) * a/mle_params(1);
                mle_params(3) = mle_params(3) * a/mle_params(1);
                mle_params(1) = a;

                [chiStat,~,~] = modelFit(this_rt,[mle_params s],nSamples,tau,max_rt);
                
                if chiStat ~= realmax
                    initial_fit = false;
                end
            
                if display_subjects
                    % Visualise initial guess
                    sim_rt = simulateProcess([mle_params s],nSamples,max_rt,tau);
                    vis_sim_newversion(this_rt,sim_rt);
                    title(['Subject ' num2str(subject) ': Initial guess (' condition_names{con} ')'])
                end
            
            end

            % Parameter estimation
            minvals(1) = a;
            maxvals(1) = a;
            loss = @(x)modelFit(this_rt,[x s],nSamples,tau,max_rt);
            nTries = 4;
            targetChiSquare = 26.1;

            cc = 0;
            carryon = true;
            chilog = nan(nTries,length(start_params)+1);
            while carryon
                [opt_params, chiSquare] = fminsearchbnd(loss,mle_params,minvals,maxvals,...
                    optimset('TolFun',0.5,'TolX',0.05,'MaxFunEvals',150,'Display','iter'));
                cc = cc + 1;
                chilog(cc,:) = [chiSquare opt_params];
                if cc >= nTries || chiSquare < targetChiSquare
                    carryon = false;
                end
                disp(['(((((((( TRY #' num2str(cc) '))))))))) ----- Subject ' num2str(subject) ', c' num2str(con) ', n ' num2str(n)])
            end
            try
                opt_params = chilog(find(chilog(:,1) == min(chilog(:,1))),2:end);
                opt_params(3) = opt_params(3)/a; % doesn't really do anything because a is 1 (X/1 = X)...
                opt_params(5) = opt_params(5)/a;

                if min(chilog(:,1)) < 80
                    
                    all_chi_square(subj,con,n) = chilog(find(chilog(:,1) == min(chilog(:,1))),1);

                    % Visualise optimal model
                    sim_rt = simulateProcess([opt_params s],nSamples,max_rt,tau);
                    if display_subjects
                        vis_sim_newversion(this_rt,sim_rt);
                        title(['Subject ' num2str(subject) ': optimised (' condition_names{con} ')'])
                        drawnow;
                    end

                    % Save
                    all_sim_data{subj,con} = [all_sim_data{subj,con}, sim_rt];
                    all_obs_data{subj,con} = [all_obs_data{subj,con}, this_rt];
                    parameter_estimates(subj,con,:) = nanmean([squeeze(parameter_estimates(subj,con,:)), opt_params'],2);
                end

            catch
                error('Chi Square too high!')
            end
            if all_chi_square(subj,con,n) > 80
                error('Chi Square too high!')
            end        
        end
        
    end
     
%     disp(['--- SUBJECT ' num2str(subj) ':'])
%     parameter_names
%     squeeze(parameter_estimates(subj,:,:))
%     
%     figure;
%     bar(squeeze(parameter_estimates(subj,:,:))')
%     set(gca,'XTickLabels',parameter_names)
%     title(['Parameter Estimates for Subject ' num2str(subj) ' (chi square = ' num2str(round(min(all_chi_square(subj,:))))...
%         ' to ' num2str(round(max(all_chi_square(subj,:)))) ')'])
%     drawnow;
    
    toc % takes about 30 minutes per subject
    
end

this_version = 2;
save(['exp' num2str(experiment) '_subjectcon_estimation_v' num2str(this_version) '.mat'],'parameter_estimates','all_chi_square','all_sim_data','all_obs_data');

%% Plot group

plot_subjects = subjects;

% Overall fit (all trials, all subjects)
Q = [.05:.05:.95];
alpha = .6;
figure
for con = 1:nCon
    
    avq_obs = [];
    avq_sim = [];
    for subj = 1:length(plot_subjects)
        avq_obs(subj,:) = quantile(mean(all_obs_data{plot_subjects(subj),con},2),Q);
        avq_sim(subj,:) = quantile(mean(all_sim_data{plot_subjects(subj),con},2),Q);
    end
        
%     for q = 1:length(Q)
%         mu = mean(avq_obs(:,q),1);
%         sigma = std(avq_obs(:,q),1)/sqrt(size(avq_obs,1));
%         plot([Q(q) Q(q)],[mu+sigma mu-sigma],'Color',colours(con,:)); hold on;
%         mu = mean(avq_sim(:,q),1);
%         sigma = std(avq_sim(:,q),1)/sqrt(size(avq_obs,1));
%         plot([Q(q) Q(q)],[mu+sigma mu-sigma],'Color',colours(con,:)); hold on;
%     end
    
    scatter(Q,mean(avq_obs,1),50,'filled','MarkerFaceColor',colours(con,:),'Marker','o','MarkerFaceAlpha',alpha,'MarkerEdgeColor','k'); hold on
    scatter(Q,mean(avq_sim,1),50,'filled','MarkerFaceColor',colours(con,:),'Marker','^','MarkerFaceAlpha',alpha,'MarkerEdgeColor','k'); hold on
    
    P = plot(Q,mean(avq_obs,1),'Color',colours(con,:),'LineWidth',1.5); hold on
%     P.Color(4) = alpha;
    P = plot(Q,mean(avq_sim,1),'Color',colours(con,:),'LineWidth',1.5,'LineStyle','--'); hold on
%     P.Color(4) = alpha;
    
end
% legend(condition_names)
title('Observed & Simulated RT Quantiles (averaged across subjects)')
ylabel('RT (seconds)')
xlabel('Quantile (%)')

% Parameter estimates (x axis = parameter)
plot_estimates = parameter_estimates(plot_subjects,:,:);
% z_estimates = zscore(plot_estimates(:,:,2:end)); % remove outliers
% for subj = 1:size(plot_estimates,1)
%     x = squeeze(z_estimates(subj,:,:));
%     if any(x(:) > 3 | x(:) < -3)
%         disp(['Subject ' num2str(plot_subjects(subj)) ' is an outlier'])
%         plot_estimates(subj,:,:) = [];
%     end
% end

figure
mu = squeeze(nanmean(plot_estimates,1));
ub = mu + squeeze(nanstd(plot_estimates,1))/sqrt(size(plot_estimates,1));
lb = mu - squeeze(nanstd(plot_estimates,1))/sqrt(size(plot_estimates,1));
jitter = linspace(-.3,.3,nCon);
for p = 1:size(plot_estimates,3)
    plot(p+jitter,mu(:,p),'Color',[.3 .3 .3],'LineWidth',1.5); hold on 
end
for subj = 1:size(plot_estimates,1)
    for p = 1:size(plot_estimates,3)
        P = plot(p+jitter,squeeze(plot_estimates(subj,:,p)),'k'); hold on
        P.Color(4) = .3;
    end
end
for con = 1:nCon
    for p = 1:size(plot_estimates,3)
        scatter(ones(size(plot_estimates,1),1)*p+jitter(con),squeeze(plot_estimates(:,con,p)),'filled','MarkerFaceColor',colours(con,:),'MarkerFaceAlpha',.5); hold on
        plot([p p]+jitter(con),[lb(con,p) ub(con,p)],'k','LineWidth',1.5); hold on 
    end
    scatter((1:length(parameter_names))+jitter(con),mu(con,:),80,'filled','MarkerFaceColor',colours(con,:),'MarkerEdgeColor','k','LineWidth',1.5); hold on
end
xlim([.5 size(plot_estimates,3)+.5])
set(gca,'XTick',1:size(plot_estimates,3))
set(gca,'XTickLabels',parameter_names)
xlabel('Parameter Type')
ylabel('Estimated Value')
title('Parameter Estimates across Group')
set(gca,'TickLength',[0 0])

% Box plots
boxdata = [];
for p = [5 2 3 4]
    boxdata = [boxdata, squeeze(plot_estimates(:,:,p))];
end

figure
boxplot(boxdata,'colors',repmat(colours,4,1),'boxstyle','filled','notch','off',...
    'symbol','k+','widths',.75,'extrememode','compress');
ylim([0 3])
set(gca,'TickLength',[0 0])
h = findobj(gcf,'tag','Outliers');
idx = flipud(repmat(1:4,1,4)');
for i = 1:length(h)
    h(i).MarkerEdgeColor = colours(idx(i,1),:);
end


spss_data = [];
for p = [2 5]
    spss_data = [spss_data, squeeze(parameter_estimates(:,:,p))];
end
