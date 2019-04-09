%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% BREAKING CONTINUOUS FLASH SUPPRESSION %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    Jessica McFadyen 2017    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

% NOTE: This script requires Psychtoolbox 3 is installed: http://psychtoolbox.org
% EEG TRIGGER CODES: 
% 1 = face onset: expected neutral
% 2 = face onset: unexpected neutral
% 3 = face onset: expected fearful
% 4 = face onset: unexpected fearful
% 5 = no mask, face onset: expected neutral
% 6 = no mask, face onset: unexpected neutral
% 7 = no mask, face onset: expected fearful
% 8 = no mask, face onset: unexpected fearful
% 9 = practice, face onset: expected neutral
% 10 = practice, face onset: unexpected neutral
% 11 = practice, face onset: expected fearful
% 12 = practice, face onset: unexpected fearful
% 99 = face offset
% responses are the same but plus 20

%% PARAMETERS ------------------------------------------------------------

screenTransparent = 0; % 0 for opaque, 1 for transparency
screen            = 2; % screen number
debugmode         = 0; % 0 for run experiment, 1 for experiment set-up only
isEEG             = 0; % 1 if there is a parallel port for EEG
dummy_tex         = 0; % 0 to use images, 1 to use black/white boxes

rand('seed',sum(100*clock)); % set up random number generator

% Subject-Specific Input
if ~debugmode
    rivalry_conf      = input('Configure rivlary distance? (1 for yes, 0 for no): ');
    subject_num       = input('Subject number?: ');
    if exist(['Output\s' num2str(subject_num) '_trial_info.mat'])
        warning('WARNING! Subject file already exists!');
    end
    subject_age       = input('Subject age?: ');
    subject_gender    = input('Subject gender? (1 for male, 2 for female): ');
    dominant_eye      = input('Dominant eye? (1 for left, 2 for right): ');
%     is_practice       = input('Do practice? (1 for yes, 0 for no) ');
    is_practice = 0;
else
    subject_num = 99;
    subject_gender = 99;
    subject_age = 99;
    dominant_eye = 1;
    block_start = randperm(2, 1);
    is_practice = 0;
end

block_start  = (mod(subject_num,2) == 0) + 1; % odd subject number = neutral first, even = fearful first

% Experiment
blocks             = 14; % no . of blocks
num_trials         = 90; % no. of trials per block
prac_blocks        = 0; % no. of practice blocks
prac_trials        = 30; % no. of trials per practice block
ITI                = 0.25:0.05:0.5; % duration of fixation cross at beginning of trial
probability        = [(5/6) (1/6)]; % highly probable, lowly probable
min_gap            = 2; % minimum no. of standards between deviants (use "pseudoRandSpace.m" plot function to check suitability)
stim_dur           = 3; % time for target to reach full contrast (in seconds)
border_width       = 1.4;
dom_size           = 1.25;
mask_control       = 0; % 0 for always show mask, 1 for set random chance of not showing mask
if mask_control
    mask_control_prob = 1; % percent of trials to hide mask
end

mondrian_contrast = 1;
min_face_contrast = 0; 
max_face_contrast = 1; 

prob_buffer        = 2;        % if block-based cueing, make sure the first N trials of each block are probable
if mod(blocks,2) == 1  % check that no. of blocks is EVEN
    warning('WARNING: Block numbers are not even! Need equal numbers Fearful and Neutral blocks.');
    go_on = input('Continue? (1 for yes, 2 for no) ');
    if go_on == 2
        sca;
        return;
    end
end

orientation_prob   = [.5 .5];  % probability of face being tilted left or right
rotation_degree    = 5;       % how many degrees the faces are rotated

% Stimuli
stim_size           = 365; % height of face stimuli
from_centre         = 293;       % distance of each centre of stimulus from centre of screen
if from_centre*2 < stim_size(1)
    error('ERROR: Left and right stimuli will overlap! Increase "from_centre" variable.');
end

desiredMonRate      = [60]; % monitor refresh rate ([] to not check)
mask_flicker        = 10;     % flicker of mondrian in Hz - 20 Hz
face_flicker        = 0;

% Display
fullscreen        = 1;         % 0 for no, 1 for yes
resolution        = [1920 1080]; 
background_colour = [.505 .505 .505];  % background colour
foreground_colour = [0 0 0];   % foreground (i.e. text) colour
fontname          = 'Calibri';
fontsize          = 40;
edge_width        = 100; % pixels from edge/centre (depending on dominant eye) where text is drawn from

if isEEG
    % Configure parallel port
    eegportaddr = hex2dec('D050'); % old PC = '378', new PC = 'D050'
    eegtriglength = 0.002; % trigger length in secs
    ioObj = io64(); % initialize the inpout32.dll system driver - change to 64 if using 64-bit computer (find "io64")
    status = io64(ioObj);
    io64(ioObj, eegportaddr, 0); % output command - set to 0/off
end

%% EXPERIMENT CONFIGURATION ---------------------------------------------

trial_info = [];
trial_info.subject.num = subject_num;
trial_info.subject.age = subject_age;
trial_info.subject.gender = subject_gender;
trial_info.dominanteye = dominant_eye;

% Load images
folder_names = {'NeutralMale','FearfulMale'; 'NeutralFemale','FearfulFemale'};
imgList = {};
for g = 1:size(folder_names,1) % rows
    for e = 1:size(folder_names,2) % cols
        temp = dir(['Stimuli/Faces/' folder_names{g,e} '/*.png']);
        for i = 1:length(temp)
            imgList{g,e}{i} = {['Stimuli/Faces/' folder_names{g,e} '/' temp(i).name]};
        end
    end
end

masks = dir('Stimuli/Masks/*.png');
masks = extractfield(masks,'name');

% Configure trial order & load images
if is_practice
    [Ptrial_order, Pimage_order, Porientation_order, Pblock_order] = configure_trials_blockbased(prac_trials,probability,min_gap,prac_blocks,imgList,orientation_prob,prob_buffer,block_start);
    trial_info.practice_trials = Ptrial_order;
    trial_info.practice_orientation = Porientation_order;
end
[trial_order, image_order, orientation_order, block_order] = configure_trials_blockbased(num_trials,probability,min_gap,blocks,imgList,orientation_prob,prob_buffer,block_start);

trial_info.trials = trial_order;
trial_info.orientation = orientation_order;

if mask_control
    trial_info.mask_control = zeros(size(trial_info.trials,1),size(trial_info.trials,2));
    for b = 1:blocks
        conditions = reshape(unique(trial_info.trials(b,:)),2,2); % first two will be standards, second two will be oddballs
        for c = 1:size(conditions,2)
           idx = [find(trial_info.trials(b,:) == conditions(1,c) | trial_info.trials(b,:) == conditions(2,c))];
           idx = datasample(idx,round(length(idx)*mask_control_prob),'replace',false);
           trial_info.mask_control(b,idx) = 1;
        end
    end
end

%% PSYCHTOOLBOX CONFIGURATION ---------------------------------------------
if ~debugmode
    sca;
    AssertOpenGL;
    PsychDefaultSetup(2); % configure OpenGL (0), keyboard (1), and screen (2)
    
    if screenTransparent == 1
        PsychDebugWindowConfiguration; % creates semi-transparent window
    end
    
    % open fullscreen window
    Screen('Preference', 'SkipSyncTests', 0);
    Screen('CloseAll')
    [window, windowRect] = PsychImaging('OpenWindow',screen,background_colour);
    Priority(MaxPriority(window)); % Switch to realtime:
    
    % check refresh rate
    monRate = Screen('NominalFrameRate', window);
    if ~isempty(desiredMonRate)
        if round(monRate)~=desiredMonRate
            error(['Change refresh rate to ' num2str(desiredMonRate) ' Hz'])
            sca;
        end
    end
    
    % check screen size
    if ~isequal(windowRect, [0 0 1920 1080])
        error('Change screen resolution to 1920 x 1080')
    end
    
    Screen('BlendFunction',window,'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % CAN SLOW DOWN COMPUTER, TRY NOT TO USE!!
    ifi = Screen('GetFlipInterval',window); % frame duration
    Screen('TextSize',window,fontsize); % set text size
    Screen('TextStyle',window,1); % 0 = normal, 1 = bold
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
    [xCentre yCentre]    = RectCenter(windowRect); % get centre coordinate
    left_centre = xCentre-from_centre;
    right_centre = xCentre+from_centre;
    
    WaitSecs('UntilTime', GetSecs+.1); % initialise WaitSecs and GetSecs
    
    % KEYBOARD
    KbName('UnifyKeyNames');
    leftKey   = KbName('LeftArrow');
    rightKey  = KbName('RightArrow');
    upKey     = KbName('UpArrow');
    downKey   = KbName('DownArrow');
    escapeKey = KbName('ESCAPE');
    spaceKey  = KbName('SPACE');
    ListenChar(0); % enables PsychHID to work (for some unknown reason)
end

% Text position
if dominant_eye == 1 % left
    text_pos = edge_width;
elseif dominant_eye == 2 % right
    text_pos = xCentre+edge_width;
end

%% RUN EXPERIMENT ---------------------------------------------------------

if ~debugmode
    HideCursor;
end

%% Load images
if ~debugmode
    [tex, mon_load] = load_stimuli(text_pos,trial_order,image_order,masks,window);
end
trial_info.images = image_order;

%% Set up stimuli distance

from_centre = 310;
vertical_pos = 547;

if ~debugmode
    if rivalry_conf
        configure_rivlary;
        disp(['X pos = ' num2str(from_centre), ', Y pos = ' num2str(vertical_pos)]);
        trial_info.vertical_distance = vertical_pos;
    end
end
trial_info.centre_distance = from_centre;
trial_info.vertical_distance = vertical_pos;

%% Trials
if ~debugmode
    
    if is_practice
        
        practice_now = 1;
        
        % Set up EEG trigger codes
        if isEEG
            eegtrig = trial_info.practice_trials;
            eegtrig(eegtrig == 31 | eegtrig == 32) = 9;
            eegtrig(eegtrig == 21 | eegtrig == 22) = 10;
            eegtrig(eegtrig == 11 | eegtrig == 12) = 11;
            eegtrig(eegtrig == 41 | eegtrig == 42) = 12;
            trial_info.practice_eegtrig = eegtrig;
        end
        
        DrawFormattedText(window,['Press any key to begin practice'],text_pos,'center');
        Screen('Flip',window);
        KbWait([],2);
        
        for b = 1:prac_blocks
            
            start_block = GetSecs;
                       
            % Run Trials
            trials;
            
            response_trials = find(~isnan(trial_info.practice_acc)); % trials where a response was made
            
            trial_info.practice_block_summary(b).response_trials = length(response_trials)/prac_trials;
            trial_info.practice_block_summary(b).acc = nanmean(trial_info.practice_acc(response_trials));
            trial_info.practice_block_summary(b).rt = nanmean(trial_info.practice_rt(response_trials));
            
            disp(['Response trials: ' num2str(length(response_trials)) ' out of ' num2str(prac_trials)]);
            disp(['Accuracy: ' num2str(mean(trial_info.practice_acc(response_trials))*100) '%']);
            disp(['RT: ' num2str(mean(trial_info.practice_rt(response_trials))) 's']);
            
            Screen('FillRect', window, background_colour);
            DrawFormattedText(window,['Accuracy = ' num2str(round(trial_info.practice_block_summary(b).acc*100)) '% \n\n End of practice block ' num2str(b) ' of ' num2str(prac_blocks)],text_pos,'center');
            Screen('Flip',window);
            
            save(['Output\s' num2str(subject_num) '_trial_info.mat'],'trial_info');
            
            end_block = 0;
            while end_block == 0
                [keyIsDown,secs,keyCode] = PsychHID('KbCheck',[],[]);
                if find(keyCode) == spaceKey
                    end_block = 1;
                elseif find(keyCode) == 67 % 'C' key
                    configure_rivalry;
                end
            end
            
            try
                end_block = GetSecs - start_block;
                trial_info.block_summary.time(b) = end_block/60; % in minutes
            end
        end
        
    end
    
    DrawFormattedText(window,['Press any key to begin experiment'],text_pos,'center');
    Screen('Flip',window);
    KbWait([],2);
    
    practice_now = 0;
    
    % Set up EEG trigger codes
    if isEEG
        eegtrig = trial_info.trials;
%         cidx = trial_info.mask_control;
        eegtrig((eegtrig == 31) | (eegtrig == 32)) = 1;
        eegtrig((eegtrig == 21) | (eegtrig == 22)) = 2;
        eegtrig((eegtrig == 11) | (eegtrig == 12)) = 3;
        eegtrig((eegtrig == 41) | (eegtrig == 42)) = 4;
%         eegtrig((eegtrig == 31 & cidx == 1) | (eegtrig == 32 & cidx == 1)) = 5;
%         eegtrig((eegtrig == 21 & cidx == 1) | (eegtrig == 22 & cidx == 1)) = 6;
%         eegtrig((eegtrig == 11 & cidx == 1) | (eegtrig == 12 & cidx == 1)) = 7;
%         eegtrig((eegtrig == 41 & cidx == 1) | (eegtrig == 42 & cidx == 1)) = 8;
        trial_info.eegtrig = eegtrig;
    end

    for b = 1:size(trial_order,1)

        start_block = GetSecs;

        % Run Trials
        trials;
        
        response_trials = find(~isnan(trial_info.acc)); % trials where a response was made

        trial_info.block_summary(b).response_trials = length(response_trials)/size(trial_order,2);
        trial_info.block_summary(b).acc = mean(trial_info.acc(response_trials));
        trial_info.block_summary(b).rt = mean(trial_info.rt(response_trials));

        disp(['Response trials: ' num2str(length(response_trials)) ' out of ' num2str(size(trial_order,2))]);
        disp(['Accuracy: ' num2str(mean(trial_info.acc(response_trials))*100) '%']);
        disp(['RT: ' num2str(mean(trial_info.rt(response_trials))) 's']);
        
        Screen('FillRect', window, background_colour);
        DrawFormattedText(window,['Accuracy = ' num2str(round(trial_info.block_summary(b).acc*100)) '% \n\n End of block ' num2str(b) ' of ' num2str(blocks)],text_pos,'center');
        Screen('Flip',window);

        save(['Output\s' num2str(subject_num) '_trial_info.mat'],'trial_info');

        end_block = 0;
        while end_block == 0
            [keyIsDown,secs,keyCode] = PsychHID('KbCheck',[],[]);
            if find(keyCode) == spaceKey
                end_block = 1;
            elseif find(keyCode) == 67 % 'C' key
                configure_rivalry;
            end
        end

        try
            end_block = GetSecs - start_block
            trial_info.block_summary.time(b) = end_block/60; % in minutes
        end
    end

    save(['Output\s' num2str(subject_num) '_workspace.mat']);

    sca
    return
end

if debugmode
    
    trial_count = trial_order(:);
    trial_count = [sum(trial_count == 31 | trial_count == 32), sum(trial_count == 41 | trial_count == 42),
        sum(trial_count == 11 | trial_count == 12), sum(trial_count == 21 | trial_count == 22)]
    
end

trials = trial_info.trials(:);
EN_idx = trials == 31 | trials == 32;
EF_idx = trials == 11 | trials == 12;
UN_idx = trials == 21 | trials == 22;
UF_idx = trials == 11 | trials == 12;

if ~debugmode
    rt = trial_info.rt(:);
    figure(1);
    rt_bar = [nanmean(rt(EN_idx)) nanmean(rt(UN_idx)); nanmean(rt(UN_idx)) nanmean(rt(UF_idx))];
    bar(rt_bar)
    title('Response Time');
    legend({'Expected','Unexpected'});

    acc = trial_info.acc(:);
    figure(2);
    acc_bar = [nanmean(acc(EN_idx)) nanmean(acc(UN_idx)); nanmean(acc(UN_idx)) nanmean(acc(UF_idx))];
    bar(acc_bar)
    title('Accuracy');
    legend({'Expected','Unexpected'});
end