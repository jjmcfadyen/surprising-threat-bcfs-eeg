%% Analyse Behavioural Data - EXP 1

close all;
clear all;
clc;

dir_data = 'D:\Scratch\bCFS_EEG_Reanalysis\data\Exp1\Behavioural'; % location of raw data files (e.g. s1_trial_info.mat)
load(fullfile(dir_data,'stai.mat')); % load state-trait anxiety scores

dir_results = 'D:\Scratch\bCFS_EEG_Reanalysis\results';

subjects            = 1:30;
remove_outliers     = 1;  % 0 for no, 1 for yes
remove_tooquick     = 1;  % 0 for no, 1 for yes
remove_tooslow      = 1;   % 0 for no, 1 for yes
tooquick            = .5; % responses shorter than 500ms are artefactual
tooslow             = 10;  % responses longer than 10s are too slow

block_include = 1:8; % blocks to analyse (i.e. exclude first 2 for learning)

% set up empty variables to be filled in for loop below
long_form = [];
OVERALL = [];
OVERALL.outliers_fast = [];
OVERALL.outliers_slow = [];
OVERALL.accuracy = [];
for s = 1:length(subjects)
    
    % Load trial info
    load(fullfile(dir_data,['s' num2str(subjects(s)) '_trial_info.mat']));
    trial_info.acc = trial_info.acc(block_include,:);
    trial_info.rt = trial_info.rt(block_include,:);
    trial_info.trials = trial_info.trials(block_include,:);
    
    % Get hits & misses
    accuracy = trial_info.acc';
    accuracy = accuracy(:);
    
    hits = ~isnan(accuracy);
    misses = isnan(accuracy);
    
    % Get responses
    try
        orientation = trial_info.orientation';
        orientation = orientation(:);
        responses = orientation;
        responses(accuracy == 0 & responses == 1) = 2;
        responses(accuracy == 0 & responses == 2) = 1;
    catch
        orientation = nan(length(trials),1);
        responses = nan(length(trials),1);
    end
    
    % Get index of condition types
    trials = trial_info.trials';
    trials = trials(:);
    
    EN = trials == 31 | trials == 32;
    UN = trials == 21 | trials == 22;
    EF = trials == 11 | trials == 12;
    UF = trials == 41 | trials == 42;             
    
    %% Calculate and plot RT (correct responses)
    
    rt = trial_info.rt';
    rt = rt(:);
    
    % Get index of how long ago the same expression was shown (i.e. how 'surprising' it is)
    surprise_idx = nan(size(trial_info.trials)); % lower number = more surprising
    for b = 1:size(trial_info.trials,1)
        
        bTrials = trial_info.trials(b,:);
        emotion_idx = bTrials == 11 | bTrials == 12 | bTrials == 41 | bTrials == 42; % 0 for neutral, 1 for fearful
        if sum(emotion_idx == 0) > sum(emotion_idx == 1) % neutral block
            deviants = 1;
        elseif sum(emotion_idx == 1) > sum(emotion_idx == 0) % fearful block
            deviants = 0;
        end
        
        % see how surprising each deviant is
        deviant_idx = find(emotion_idx == deviants);
        surprise_idx(b,deviant_idx) = [deviant_idx(1) diff(deviant_idx)];
        
        % see how many repetitions there have been of each standard
        for i = 1:length(deviant_idx)
            if i == 1
                these_trials = 1:deviant_idx(i)-1;
                rep_idx = fliplr([1:length(these_trials)] - length(these_trials));
            else
                these_trials = deviant_idx(i-1)+1:deviant_idx(i)-1;
                rep_idx = fliplr([1:length(these_trials)]- length(these_trials));
            end
            surprise_idx(b,these_trials) = rep_idx;
        end
        these_trials = deviant_idx(end)+1:size(surprise_idx,2);
        surprise_idx(b,these_trials) = fliplr([1:length(these_trials)]- length(these_trials));
        
    end
    surprise_idx = surprise_idx';
    surprise_idx = surprise_idx(:);
    surprise_idx(isnan(surprise_idx)) = 99; % so that I can find them in R later
           
    %% Save individual trial data
    
    ntrials = length(trials);
    
    ntrials_perblock = size(trial_info.trials,2);
    nblocks = size(trial_info.trials,1);
    
    nblocktrials = repmat(1:ntrials_perblock,1,nblocks)';
    nblockcount = repmat(1:nblocks,ntrials_perblock,1);
    nblockcount = nblockcount(:);
    
    emotion_idx = (EF | UF) + 1; % emotion (1 = neutral, 2 = fearful)
    expectation_idx = (UN | UF) + 1; % expectation (1 = expected, 2 = unexpected)
    
    if trial_info.trials(1,1) == 31 || trial_info.trials(1,1) == 32 || trial_info.trials(1,1) == 41 || trial_info.trials(1,1) == 42
        block_order = 1; % neutral block first
    else block_order = 2; % fearful block first
    end
    
    subject_long_form = [
                    ones(ntrials,1)*s,... % subject number
                    ones(ntrials,1)*trial_info.subject.gender,... % gender
                    ones(ntrials,1)*trial_info.subject.age,... % age
                    ones(ntrials,1)*stai.full(s,1),... % STAI additive score
                    ones(ntrials,1)*stai.state(s,1),... % STAI state score
                    ones(ntrials,1)*stai.trait(s,1),... % STAI trait score
                    ones(ntrials,1)*block_order,... % order
                    nblockcount,... % block number
                    nblocktrials,... % trial in block
                    [1:ntrials]',... % trial in exp
                    surprise_idx,... % how frequently the emotion has been presented recently (lower number = more surprising)
                    emotion_idx,... % emotion (1 = neutral, 2 = fearful)
                    expectation_idx,... % expectation (1 = expected, 2 = unexpected)
                    rt,... % rt
                    responses,... % response type (left or right)
                    orientation]; % left or right

    % Select only correct responses
    subject_long_form = subject_long_form(accuracy == 1,:);
    ntrials = size(subject_long_form,1);
    
    % Identify fast outliers
    if remove_tooquick
         outliers_fast = subject_long_form(:,14) < tooquick;
    else outliers_fast = zeros(ntrials,1);
    end
    OVERALL.outliers_fast = [OVERALL.outliers_fast; sum(outliers_fast)];
    
    % Identify extremely slow responses
    if remove_tooslow
         outliers_slow = subject_long_form(:,14) > tooslow;
    else outliers_slow = zeros(ntrials,1);
    end
    OVERALL.outliers_slow = [OVERALL.outliers_slow; sum(outliers_slow)];
               
    outliers = outliers_slow | outliers_fast;
    
    % Remove the outliers
    if remove_outliers
        subject_long_form = subject_long_form(outliers == 0,:);
    end

    % Pool together
    long_form = [long_form; subject_long_form];
    
    % Accuracy
    OVERALL.accuracy(s,1) = mean(accuracy(EN & (rt > tooquick & rt < tooslow)));
    OVERALL.accuracy(s,2) = mean(accuracy(UN & (rt > tooquick & rt < tooslow)));
    OVERALL.accuracy(s,3) = mean(accuracy(EF & (rt > tooquick & rt < tooslow)));
    OVERALL.accuracy(s,4) = mean(accuracy(UF & (rt > tooquick & rt < tooslow)));
    
end

% Save trial_data (to be used in LME analysis in R)
header_names = {'Subject','Gender','Age','STAI','State','Trait',...
                'Order','Block','BlockTrial','ExpTrial','SurpriseIdx','Emotion','Expectation','RT','Response','Orientation'};        
            
commaHeader = [header_names; repmat({','}, 1, numel(header_names))]; %insert commaas
textHeader = cell2mat(commaHeader(:)'); %cHeader in text with commas

filename = fullfile(dir_results,'trial_data_exp1.csv');

fid = fopen(filename,'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);

dlmwrite(filename,long_form,'-append');

% Get trial counts
emotion_col = 12;
expectation_col = 13;

cidx = {};
cidx{1} = long_form(:,emotion_col) == 1 & long_form(:,expectation_col) == 1; % EN
cidx{2} = long_form(:,emotion_col) == 1 & long_form(:,expectation_col) == 2; % UN
cidx{3} = long_form(:,emotion_col) == 2 & long_form(:,expectation_col) == 1; % EF
cidx{4} = long_form(:,emotion_col) == 2 & long_form(:,expectation_col) == 2; % UF

trial_counts = [];
for s = 1:length(subjects)
    for c = 1:4
        trial_counts(s,c) = sum(long_form(:,1) == s & cidx{c});
    end
end
