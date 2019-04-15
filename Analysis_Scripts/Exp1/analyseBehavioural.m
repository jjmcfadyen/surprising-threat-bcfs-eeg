%% Analyse Behavioural Data - EXP 1

close all;
clear all;
clc;

dir_data = 'Data\Exp1\Behavioural'; % location of raw data files (e.g. s1_trial_info.mat)
load(fullfile(dir_data,'stai.mat')); % load state-trait anxiety scores

subjects            = 1:30;
only_after_dev      = 0; % 0 for no, X no. of standards afer deviant
dev_pos             = 0; % 0 to ignore, 1 to code for deviant position (i.e. how many standards since last deviant)
remove_outliers     = 1; % 0 for no, 1 for yes
remove_tooquick     = 1; % 0 for no, 1 for yes
tooquick            = .5; % responses shorter than 200ms are artefactual
outlier_def         = 5; % no. of std away from mean
norm_check          = 0; % 1 for yes, 0 for no - uses Kolmogorov-Smirnov test
error_bars          = 1; % 0 for none, 1 for standard deviation, 2 for standard error
av_type             = 1; % 1 for mean, 2 for median

block_include = 3:8; % blocks to analyse (i.e. exclude first 2 for learning)
plot_subjects = 0;

% set up empty variables to be filled in for loop below
group_data = [];
group_trial_count = [];
group_dev_pos.UN = NaN(length(subjects),101);
group_dev_pos.UF = NaN(length(subjects),101);
group_outliers = [];
group_outliers_fast = [];
trial_data = [];
for s = 1:length(subjects)
    
    % Load trial info
    load(fullfile(dir_data,['s' num2str(subjects(s)) '_trial_info.mat']));
    trial_info.acc = trial_info.acc(block_include,:);
    trial_info.rt = trial_info.rt(block_include,:);
    trial_info.trials = trial_info.trials(block_include,:);
    
    % Get trial counts (~isnan)
    hits = ~isnan(trial_info.acc(:));
    misses = isnan(trial_info.acc(:));
    
    % Get index of condition types
    trials = trial_info.trials';
    trials = trials(:);
    
    EN = trials == 31 | trials == 32;
    UN = trials == 21 | trials == 22;
    EF = trials == 11 | trials == 12;
    UF = trials == 41 | trials == 42;
    
    trial_count = [sum(EN(hits)), sum(EF(hits));
                   sum(UN(hits)), sum(UF(hits))];               
    nano_count = [sum(EN(misses)), sum(EF(misses));
                   sum(UN(misses)), sum(UF(misses))];
    
    accuracy = trial_info.acc';
    accuracy = accuracy(:);
    
    %% Calculate and plot RT (correct responses)
    
    rt = trial_info.rt';
    rt = rt(:);
    correct_responses = accuracy == 1;
    
    % Recode standards
    recoded_trials = trials;
    if only_after_dev > 0
        
        EN = zeros(length(trials),1);
        idx = find(UF)+only_after_dev;
        EN(idx,1) = true;
        if length(EN) > length(recoded_trials)
            EN = EN(1:length(recoded_trials),1);
        end
        EN = EN & ~UF & ~EF & ~UN;
        
        EF = zeros(length(trials),1);
        idx = find(UN)+only_after_dev;
        EF(idx,1) = true;
        if length(EF) > length(recoded_trials)
            EF = EF(1:length(recoded_trials),1);
        end
        EF = EF & ~EN & ~UN & ~UF;
        
        trial_count = [sum(EN(hits)), sum(EF(hits));...
                   sum(UN(hits)), sum(UF(hits))];
        group_trial_count(s,:) = [sum(EN(hits)) sum(UN(hits)) sum(EF(hits)) sum(UF(hits))];
    end
       
    % Remove outliers
    if remove_tooquick
         outliers_fast = rt < tooquick;
    else outliers_fast = zeros(length(rt),1);
    end
    accuracy = accuracy(outliers_fast == 0,1);
    group_outliers_fast = [group_outliers_fast; sum(outliers_fast)];
    
    outliers = rt > (mean(rt(correct_responses)) + outlier_def*(std(rt(correct_responses)))); % | ...
               %rt < (mean(rt(correct_responses)) - outlier_def*(std(rt(correct_responses))));
    group_outliers = [group_outliers; sum(outliers)];
               
    outliers = outliers | outliers_fast;
           
    if ~remove_outliers
        scatter_EN = rt(EN & correct_responses);
        scatter_EF = rt(EF & correct_responses);
        scatter_UN = rt(UN & correct_responses);
        scatter_UF = rt(UF & correct_responses);
        
        accuracy_bar = [nanmean(accuracy(EN)), nanmean(accuracy(EF));
                        nanmean(accuracy(UN)), nanmean(accuracy(UF))];
        
    elseif remove_outliers
        
        scatter_EN = rt(EN & correct_responses & ~outliers);
        scatter_EF = rt(EF & correct_responses & ~outliers);
        scatter_UN = rt(UN & correct_responses & ~outliers);
        scatter_UF = rt(UF & correct_responses & ~outliers);
        
        accuracy_bar = [nanmean(accuracy(EN(~outliers))), nanmean(accuracy(EF(~outliers)));
                        nanmean(accuracy(UN(~outliers))), nanmean(accuracy(UF(~outliers)))];
        
    end
     

    if plot_subjects
        figure;
        
        B1 = bar(accuracy_bar);
        set(B1(1),'FaceColor','b');
        set(B1(2),'FaceColor','r');
        legend('Neutral','Fearful');
        set(gca,'XTickLabel',{'Expected','Unexpected'});
        ylabel('Accuracy');

        saveas(gcf,['figures\accuracy_s' num2str(s)],'jpeg');
    end
    
    %% Check for normality
    if norm_check
        
        if remove_outliers
            rt_data = rt(correct_responses & ~outliers);
        else rt_data = rt(correct_responses);
        end
        
        subplot(1,3,1);
        histfit(rt_data);
        [H,P] = kstest(rt_data);
        title(['Untransformed: H = ' num2str(H)]);
        
        transformed_rt = log10(rt_data);
        subplot(1,3,2);
        histfit(transformed_rt);
        [H,P] = kstest(transformed_rt);
        title(['Log Transformed: H = ' num2str(H)]);
        
        if ~remove_outliers
            scatter_EN = transformed_rt(EN(correct_responses));
            scatter_EF = transformed_rt(EF(correct_responses));
            scatter_UN = transformed_rt(UN(correct_responses));
            scatter_UF = transformed_rt(UF(correct_responses));
        else
            scatter_EN = transformed_rt(EN(correct_responses & ~outliers));
            scatter_EF = transformed_rt(EF(correct_responses & ~outliers));
            scatter_UN = transformed_rt(UN(correct_responses & ~outliers));
            scatter_UF = transformed_rt(UF(correct_responses & ~outliers));
        end
        
    end
    
    if av_type == 1
        rt_bar = [mean(scatter_EN), mean(scatter_UN);
                    mean(scatter_EF), mean(scatter_UF)];
    elseif av_type == 2
        rt_bar = [median(scatter_EN), median(scatter_UN);
                    median(scatter_EF), median(scatter_UF)];
    end
    
    if norm_check  & plot_subjects    
        subplot(1,3,3);
    end
    
    if plot_subjects
        % Plot bar graph
        B2 = bar(rt_bar);
        set(B2(1),'FaceColor','w');
        set(B2(2),'FaceColor','w');
        set(B2(1),'EdgeColor','b');
        set(B2(1),'LineWidth',2);
        set(B2(2),'EdgeColor','r');
        set(B2(2),'LineWidth',2);
        hold on;

        % Plot error bars
        if error_bars == 1
            upper_lim = [mean(scatter_EN)+(std(scatter_EN)/2), mean(scatter_UN)+(std(scatter_UN)/2),...
                         mean(scatter_EF)+(std(scatter_EF)/2), mean(scatter_UF)+(std(scatter_UF)/2)];
            lower_lim = [mean(scatter_EN)-(std(scatter_EN)/2), mean(scatter_UN)-(std(scatter_UN)/2),...
                         mean(scatter_EF)-(std(scatter_EF)/2), mean(scatter_UF)-(std(scatter_UF)/2)];
        elseif error_bars == 2
            upper_lim = [mean(scatter_EN)+((std(scatter_EN)/sqrt(length(scatter_EN)))/2), mean(scatter_UN)+((std(scatter_UN)/sqrt(length(scatter_UN)))/2),...
                         mean(scatter_EF)+((std(scatter_EF)/sqrt(length(scatter_EF)))/2), mean(scatter_UF)+((std(scatter_UF)/sqrt(length(scatter_EN)))/2)];
            lower_lim = [mean(scatter_EN)-((std(scatter_EN)/sqrt(length(scatter_EN)))/2), mean(scatter_UN)-((std(scatter_UN)/sqrt(length(scatter_UN)))/2),...
                         mean(scatter_EF)-((std(scatter_EF)/sqrt(length(scatter_EF)))/2), mean(scatter_UF)-((std(scatter_UF)/sqrt(length(scatter_EN)))/2)];
        end

        errorbaridx = [.8 .9 .85; 1.1 1.2 1.15; 1.8 1.9 1.85; 2.1 2.2 2.15];   
        for i = 1:size(errorbaridx,1)
            plot([errorbaridx(i,1) errorbaridx(i,2)],repmat(upper_lim(i),1,2),'k'); hold on;
            plot([errorbaridx(i,1) errorbaridx(i,2)],repmat(lower_lim(i),1,2),'k'); hold on;
            plot([errorbaridx(i,3) errorbaridx(i,3)],[upper_lim(i) lower_lim(i)],'k'); hold on;
        end

        % Plot actual trials
        dot_size = 20;
        S1 = scatter(ones(1,length(scatter_EN))*.85,scatter_EN,dot_size,'b','Filled'); hold on;
        S2 = scatter(ones(1,length(scatter_UN))*1.15,scatter_UN,dot_size,'r','Filled'); hold on;
        S3 = scatter(ones(1,length(scatter_EF))*1.85,scatter_EF,dot_size,'b','Filled'); hold on;
        S4 = scatter(ones(1,length(scatter_UF))*2.15,scatter_UF,dot_size,'r','Filled');

        legend('Expected','Unexpected');
        set(gca,'XTickLabel',{'Neutral','Fearful'});
        title('RT (correct responses)');

        saveas(gcf,['figures\rt_c' num2str(s)],'jpeg');
    end
    
    % Find position of deviants
    if dev_pos
        
        UN_pos = find(UN);
        UN_pos = [UN_pos, [UN_pos(1,1); diff(UN_pos)]];
        UF_pos = find(UF);
        UF_pos = [UF_pos, [UF_pos(1,1); diff(UF_pos)]];
        
        unique_UN_pos = unique(UN_pos(:,2));
        for i = 1:length(unique_UN_pos)
            group_dev_pos.UN(s,unique_UN_pos(i)) = mean(rt(find(UN_pos(:,2) == unique_UN_pos(i) & correct_responses(UN) & ~outliers(UN))));
        end
        
        unique_UF_pos = unique(UF_pos(:,2));
        for i = 1:length(unique_UF_pos)
            group_dev_pos.UF(s,unique_UF_pos(i)) = mean(rt(find(UF_pos(:,2) == unique_UF_pos(i) & correct_responses(UF) & ~outliers(UF))));
        end
        
    end
    
    %% Calculate significance of RT
    
%     disp('-------------------------------');
%     disp(['        RT FOR SUBJECT ' num2str(subjects(s))]);
%     disp('-------------------------------');
%     disp('--- EMOTION ---');
%     [~,p] = ttest2(scatter_EN,scatter_EF);
%     disp(['Expected Neutral vs. Expected Fearful: p = ' num2str(p)]);
%     
%     [~,p] = ttest2(scatter_UN,scatter_UF);
%     disp(['Unexpected Neutral vs. Unexpected Fearful: p = ' num2str(p)]);
%     
%     disp('--- PREDICTION ---');
%     [~,p] = ttest2(scatter_EN,scatter_UN);
%     disp(['Expected Neutral vs. Unexpected Neutral: p = ' num2str(p)]);
%     
%     [~,p] = ttest2(scatter_EF,scatter_UF);
%     disp(['Expected Fearful vs. Unexpected Fearful: p = ' num2str(p)]);
%     
%     close all;   
%     

    %% Save Group Data
    
    % Demographics
    group_data(s,1) = s;
    group_data(s,2) = trial_info.subject.age;
    group_data(s,3) = trial_info.subject.gender;
    group_data(s,4) = trial_info.dominanteye;
    
    % ACC
    group_data(s,5) = accuracy_bar(1,1);
    group_data(s,6) = accuracy_bar(1,2);
    group_data(s,7) = accuracy_bar(2,1);
    group_data(s,8) = accuracy_bar(2,2);
    
    % RT
    group_data(s,9) = rt_bar(1,1);
    group_data(s,10) = rt_bar(1,2);
    group_data(s,11) = rt_bar(2,1);
    group_data(s,12) = rt_bar(2,2);
    
    group_trial_count(s,:) = [length(scatter_EN), length(scatter_UN), length(scatter_EF), length(scatter_UF)];
    
    %% Save individual trial data
    
    if ~remove_outliers
        outliers = zeros(720,1);
    end
    
    nblocktrials = repmat(1:90,1,8)';
    nblockcount = repmat(1:8,90,1);
    nblockcount = nblockcount(:);
    nblockcount = nblockcount(~outliers & correct_responses,1);
    
    ntrials = sum(~outliers & correct_responses);
    ntrials_idx = [1:(8*90)]';
    nblocktrials = nblocktrials(~outliers & correct_responses,1);

    emotion_idx = (EF | UF) + 1; % emotion (1 = neutral, 2 = fearful)
    expectation_idx = (UN | UF) + 1; % expectation (1 = expected, 2 = unexpected)
    
    if trial_info.trials(1,1) == 31 || trial_info.trials(1,1) == 32 || trial_info.trials(1,1) == 41 || trial_info.trials(1,1) == 42
        block_order = 1; % neutral block first
    else block_order = 2; % fearful block first
    end
    
    trial_data = [trial_data; ...
        
                  ones(ntrials,1)*s,... % subject number
                  ones(ntrials,1)*trial_info.subject.gender,... % gender
                  ones(ntrials,1)*trial_info.subject.age,... % age
                  ones(ntrials,1)*stai.full(s,1),... % STAI additive score
                  ones(ntrials,1)*stai.state(s,1),... % STAI state score
                  ones(ntrials,1)*stai.trait(s,1),... % STAI trait score                  
                  ones(ntrials,1)*stai.meancentred_full(s,1),... % STAI mean-centred additive score
                  ones(ntrials,1)*stai.meancentred_state(s,1),... % STAI mean-centred state score
                  ones(ntrials,1)*stai.meancentred_trait(s,1),... % STAI mean-centred trait score                
                  ones(ntrials,1)*block_order,... % order
                  nblockcount,... % block number
                  nblocktrials,... % trial in block
                  ntrials_idx(~outliers & correct_responses),... % trial in exp                  
                  emotion_idx(~outliers & correct_responses),... % emotion (1 = neutral, 2 = fearful)
                  expectation_idx(~outliers & correct_responses),... % expectation (1 = expected, 2 = unexpected)                 
                  rt(correct_responses & ~outliers)]; % rt

end

% Save trial_data (to be used in LME analysis in R)
header_names = {'Subject','Gender','Age','STAI','State','Trait','mcSTAI','mcState','mcTrait',...
                'Order','Block','ExpTrial','Emotion','Expectation','RT'};        
            
commaHeader = [header_names; repmat({','}, 1, numel(header_names))]; %insert commaas
textHeader = cell2mat(commaHeader(:)'); %cHeader in text with commas

fid = fopen('trial_data.csv','w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);

dlmwrite('trial_data.csv',trial_data,'-append');

% Get trial counts
emotion_col = 14;
expectation_col = 15;

cidx = {};
cidx{1} = trial_data(:,emotion_col) == 1 & trial_data(:,expectation_col) == 1; % EN
cidx{2} = trial_data(:,emotion_col) == 1 & trial_data(:,expectation_col) == 2; % UN
cidx{3} = trial_data(:,emotion_col) == 2 & trial_data(:,expectation_col) == 1; % EF
cidx{4} = trial_data(:,emotion_col) == 2 & trial_data(:,expectation_col) == 2; % UF

trial_counts = [];
for s = 1:30
    for c = 1:4
        trial_counts(s,c) = sum(trial_data(:,1) == s & cidx{c});
    end
end
