clear all
close all
clc

addpath('LBA-master');

%% Variables

experiment = 1;

% maximum trial length
if experiment == 1
    max_rt = 10;
elseif experiment == 2
    max_rt = 3;
end

min_rt = .5;

%% Load data

T = readtable(['D:\Scratch\bCFS_EEG_Reanalysis\results\trial_data_exp' num2str(experiment) '.csv'],'FileType','spreadsheet','ReadVariableNames',true);
T = T(T.RT > .5,:);

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

parameter_names = {'v' % mean drift rate
                   'A' % maximum start point
                   'b' % response threshold
                   'sv' % drift rate variability
                   't0'}; % nondecision time
          
%% Define Models

all_models = {};

model = [];
model.v = 1;
model.A = 1;
model.b = 1;
model.sv = 1;
model.t0 = 1;
all_models{length(all_models)+1} = model;

model = [];
model.v = nCon;
model.A = 1;
model.b = 1;
model.sv = 1;
model.t0 = 1;
all_models{length(all_models)+1} = model;

model = [];
model.v = 1;
model.A = nCon;
model.b = 1;
model.sv = 1;
model.t0 = 1;
all_models{length(all_models)+1} = model;

model = [];
model.v = 1;
model.A = 1;
model.b = nCon;
model.sv = 1;
model.t0 = 1;
all_models{length(all_models)+1} = model;

model = [];
model.v = nCon;
model.A = nCon;
model.b = 1;
model.sv = 1;
model.t0 = 1;
all_models{length(all_models)+1} = model;

model = [];
model.v = 1;
model.A = nCon;
model.b = nCon;
model.sv = 1;
model.t0 = 1;
all_models{length(all_models)+1} = model;

model = [];
model.v = nCon;
model.A = 1;
model.b = nCon;
model.sv = 1;
model.t0 = 1;
all_models{length(all_models)+1} = model;

model = [];
model.v = nCon;
model.A = nCon;
model.b = nCon;
model.sv = 1;
model.t0 = 1;
all_models{length(all_models)+1} = model;

% %% MLE parameter estimation for GROUP
% % Continuous maximum likelihood estimation of LBA model
% 
% rt = T.RT;
% data = [];
% data.rt = rt*1000;
% data.cond = conditions;
% data.response = ones(length(rt),1);
% data.stim = data.response;
% 
% data = LBA_clean(data);
% 
% all_params = nan(length(all_models),nCon,5);
% all_ll = nan(length(all_models),1);
% for m = 1:length(all_models)
%     
%     disp('---------------------------------------------------------------------')
%     disp(['Running model ' num2str(m) ' of ' num2str(length(all_models)) '...'])
%     disp('---------------------------------------------------------------------')
%     
%     model = all_models{m};
% 
%     pArray = [ones(1,model.v)*.8 ones(1,model.A)*300 ones(1,model.b)*150 ones(1,model.sv)*0.4 ones(1,model.t0)*200];
% 
%     [params ll] = LBA_mle(data, model, pArray);
% 
%     cor = data.response == data.stim;
%     [v a b sv t0] = LBA_parse(model, params, nCon);
%     
%     all_params(m,:,:) = [v a b sv t0];
%     all_ll(m,1) = ll;
%     
% end
% 
% % Plot group model fit
% figure
% bar(all_ll)
% title('Log Likelihoods')
% xlabel('Models')
% 
% % Plot group parameters
% winning_model = find(all_ll == max(all_ll));
% model = all_models{winning_model};
% figure
% for p = 1:5
%     subplot(2,3,p)
%     B = bar(all_params(winning_model,:,p));
%     set(gca,'XTickLabels',condition_names)
%     title(parameter_names{p})
% end

%% MLE parameter estimation for EACH SUBJECT
% Continuous maximum likelihood estimation of LBA model

all_params = nan(length(subjects),length(all_models),nCon,5);
all_ll = nan(length(subjects),length(all_models),1);
for subj = 1:length(subjects)
    
    rt = T.RT(T.Subject == subjects(subj));
    data = [];
    data.rt = rt*1000;
    data.cond = conditions(T.Subject == subjects(subj));
    data.response = ones(length(rt),1);
    data.stim = data.response;
    
    data = LBA_clean(data);
    
    for m = 1:length(all_models)
        
        disp('---------------------------------------------------------------------')
        disp(['--- Running subject ' num2str(subjects(subj)) ', model ' num2str(m) ' of ' num2str(length(all_models)) '...'])
        disp('---------------------------------------------------------------------')
        
        model = all_models{m};
        
        pArray = [ones(1,model.v)*.8 ones(1,model.A)*300 ones(1,model.b)*150 ones(1,model.sv)*0.4 ones(1,model.t0)*200];
        
        [params ll] = LBA_mle(data, model, pArray);
        
        cor = data.response == data.stim;
        [v a b sv t0] = LBA_parse(model, params, nCon);
        
        all_params(subj,m,:,:) = [v a b sv t0];
        all_ll(subj,m,1) = ll;
        
    end
end

% Plot model fit
figure
bar(-all_ll)
title('Log Likelihoods per Subject')
xlabel('Models')

winning_model = nan(length(subjects),2);
for subj = 1:length(subjects)
    [tmp, sort_idx] = sort(all_ll(subj,:));
    winning_model(subj,1) = sort_idx(1);
    winning_model(subj,2) = tmp(2)-tmp(1);
end
figure
histogram(winning_model(:,1))

% % Plot parameters
% winning_model = find(all_ll == max(all_ll));
% model = all_models{winning_model};
% figure
% for p = 1:5
%     subplot(2,3,p)
%     B = bar(all_params(winning_model,:,p));
%     set(gca,'XTickLabels',condition_names)
%     title(parameter_names{p})
% end
