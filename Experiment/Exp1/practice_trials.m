%% Run Trials

%% FIXATION CROSS

% Fixation cross
leftRect = CenterRectOnPointd([0 0 screenXpixels/2 screenYpixels],left_centre, vertical_pos);
rightRect = CenterRectOnPointd([0 0 screenXpixels/2 screenYpixels],right_centre, vertical_pos);
centre_line = CenterRectOnPointd([0 0 3 screenYpixels],xCentre,vertical_pos);
Screen('FillRect', window, background_colour, leftRect);
Screen('FillRect', window, background_colour, rightRect);
DrawFormattedText(window,'+','center','center',[],[],[],[],[],[],leftRect);
DrawFormattedText(window,'+','center','center',[],[],[],[],[],[],rightRect);
Screen('FillRect', window, [0 0 0], centre_line);
Screen('Flip',window);
WaitSecs(ITI(randperm(length(ITI),1)));

%% TARGET

% Which orientation?
if Porientation_order(b,t) == 1
    rotation =  -rotation_degree; % left
elseif Porientation_order(b,t) == 2
    rotation = rotation_degree; % right
end

% Draw centre line
centre_line = CenterRectOnPointd([0 0 3 screenYpixels],xCentre,vertical_pos);
Screen('FillRect', window, [0 0 0], centre_line);
Screen('Flip',window,[],1);

% Fade in target using sine wave function
if dominant_eye == 1
    domRect_inner = CenterRectOnPointd([0 0 stim_size(1)*1.25 stim_size(2)*1.25],left_centre, vertical_pos);
    subRect_inner_face = CenterRectOnPointd([0 0 stim_size(1) stim_size(2)],right_centre, vertical_pos);
    subRect_inner_grey = CenterRectOnPointd([0 0 stim_size(1)*1.25 stim_size(2)*1.25],right_centre, vertical_pos);
    domRect_outer = CenterRectOnPointd([0 0 stim_size(1)*border_width stim_size(2)*border_width],left_centre, vertical_pos);
    subRect_outer = CenterRectOnPointd([0 0 stim_size(1)*border_width stim_size(2)*border_width],right_centre, vertical_pos);
elseif dominant_eye == 2
    domRect_inner = CenterRectOnPointd([0 0 stim_size(1)*1.25 stim_size(2)*1.25],right_centre, vertical_pos);
    subRect_inner_face = CenterRectOnPointd([0 0 stim_size(1) stim_size(2)],left_centre, vertical_pos);
    subRect_inner_grey = CenterRectOnPointd([0 0 stim_size(1)*1.25 stim_size(2)*1.25],left_centre, vertical_pos);
    domRect_outer = CenterRectOnPointd([0 0 stim_size(1)*border_width stim_size(2)*border_width],right_centre, vertical_pos);
    subRect_outer = CenterRectOnPointd([0 0 stim_size(1)*border_width stim_size(2)*border_width],left_centre, vertical_pos);
end

% Present stimuli
mondrian_counter = 1;
curr_mondrian = 1;
response = 0;
frame_counter = 0;
num_frames = stim_dur/ifi;
inc = 1/num_frames;
start_trial = GetSecs; % start reaction time timer
idx = start_trial:(1/flicker_rate):(start_trial+200);
if face_flicker > 0
    fidx = [start_trial:(1/face_flicker/2):(start_trial+200)]';
    fidx(1:length(fidx),2) = 1;
    fidx(1:2:end,2) = 0;
    face_flicker_counter = 1;
end
if isEEG
    % Send trigger for trial start
    trigger_code = 100 + trial_info.trials(b,t); % first number: 1 for face, 2 for response, second & third numbers: face type
    io64(ioObj, eegportaddr, trigger_code);  % Turns trigger ON
    wait(eegtriglength);
    io64(ioObj, eegportaddr, 0);     % turns trigger OFF
end
if end_trial
    while ~response
        
        % contrast values
        face_contrast = min(max_face_contrast,min_face_contrast+(inc*frame_counter));
        if ramp_mask
            mondrian_contrast = 1-face_contrast;
        else mondrian_contrast = 1;
        end
        frame_counter = frame_counter+1;
        
        % select mondrian
        if GetSecs > idx(mondrian_counter)
            mondrian_counter = mondrian_counter + 1;
            mondrian_options = find(1:length(mondrians)~=curr_mondrian);
            curr_mondrian = mondrian_options(randi(length(mondrian_options)));
        end
        
        % draw face
        if ~face_flicker
            Screen('FillRect', window, [0 0 0], subRect_outer);
            Screen('FillRect', window, background_colour, subRect_inner_grey);
            Screen('DrawTexture',window,tex(b,t),[],subRect_inner_face,rotation,[],face_contrast);
        elseif face_flicker > 0 && GetSecs > fidx(face_flicker_counter)
            face_flicker_counter = face_flicker_counter + 1;
            if fidx(face_flicker_counter,2) == 1
                Screen('FillRect', window, [0 0 0], subRect_outer);
                Screen('FillRect', window, background_colour, subRect_inner_grey);
                Screen('DrawTexture',window,tex(b,t),[],subRect_inner_face,rotation,[],face_contrast);
            else
                Screen('FillRect', window, [0 0 0], subRect_outer);
                Screen('FillRect', window, background_colour, subRect_inner_grey);
            end
        end
        
        % draw mondrian
        Screen('FillRect', window, [0 0 0], domRect_outer);
        Screen('FillRect', window, background_colour, domRect_inner);
        Screen('DrawTexture',window,mon_load( curr_mondrian),[],domRect_inner,0,[],mondrian_contrast);
        
        % draw fixation crosses
        DrawFormattedText(window,'+','center','center',0,[],[],[],[],[],leftRect);
        DrawFormattedText(window,'+','center','center',0,[],[],[],[],[],rightRect);
        
        % flip
        Screen('Flip',window,[],1);
        
        % check for response
        [keyIsDown,secs,keyCode] = KbCheck;
        if keyCode(escapeKey)
            sca;
            return;
        elseif keyCode(leftKey)
            if isEEG
                % Send trigger for trial end
                trigger_code = 200 + trial_info.trials(b,t); % first number: 1 for face, 2 for response, second & third numbers: face type
                io64(ioObj, eegportaddr, trigger_code);  % Turns trigger ON
                wait(eegtriglength);
                io64(ioObj, eegportaddr, 0);     % turns trigger OFF
            end
            response = 1;
            response_time = GetSecs - start_trial;
        elseif keyCode(rightKey)
            if isEEG
                % Send trigger for trial end
                trigger_code = 200 + trial_info.trials(b,t); % first number: 1 for face, 2 for response, second & third numbers: face type
                io64(ioObj, eegportaddr, trigger_code);  % Turns trigger ON
                wait(eegtriglength);
                io64(ioObj, eegportaddr, 0);     % turns trigger OFF
            end
            response = 2;
            response_time = GetSecs - start_trial;
        end
        
    end
else
    while GetSecs < start_trial+stim_dur
        
        % contrast values
        face_contrast = min(max_face_contrast,min_face_contrast+(inc*frame_counter));
        if ramp_mask
            mondrian_contrast = 1-face_contrast;
        else mondrian_contrast = 1;
        end
        frame_counter = frame_counter+1;
        
        % select mondrian
        if GetSecs > idx(mondrian_counter)
            mondrian_counter = mondrian_counter + 1;
            mondrian_options = find(1:length(mondrians)~=curr_mondrian);
            curr_mondrian = mondrian_options(randi(length(mondrian_options)));
        end
        
        % draw face
        if ~face_flicker
            Screen('FillRect', window, [0 0 0], subRect_outer);
            Screen('FillRect', window, background_colour, subRect_inner_grey);
            Screen('DrawTexture',window,tex(b,t),[],subRect_inner_face,rotation,[],face_contrast);
        elseif face_flicker > 0 && GetSecs > fidx(face_flicker_counter)
            face_flicker_counter = face_flicker_counter + 1;
            if fidx(face_flicker_counter,2) == 1
                Screen('FillRect', window, [0 0 0], subRect_outer);
                Screen('FillRect', window, background_colour, subRect_inner_grey);
                Screen('DrawTexture',window,tex(b,t),[],subRect_inner_face,rotation,[],face_contrast);
            else
                Screen('FillRect', window, [0 0 0], subRect_outer);
                Screen('FillRect', window, background_colour, subRect_inner_grey);
            end
        end
        
        % draw mondrian
        Screen('FillRect', window, [0 0 0], domRect_outer);
        Screen('FillRect', window, background_colour, domRect_inner);
        Screen('DrawTexture',window,mon_load( curr_mondrian),[],domRect_inner,0,[],mondrian_contrast);
        
        % draw fixation crosses
        DrawFormattedText(window,'+','center','center',0,[],[],[],[],[],leftRect);
        DrawFormattedText(window,'+','center','center',0,[],[],[],[],[],rightRect);
        
        % flip
        Screen('Flip',window,[],1);
        
        % check for response
        [keyIsDown,secs,keyCode] = KbCheck;
        if keyCode(escapeKey)
            sca;
            return;
        elseif keyCode(leftKey)
            if isEEG
                % Send trigger for trial end
                trigger_code = 400 + trial_info.trials(b,t); % first number: 1 for face, 2 for response, second & third numbers: face type
                io64(ioObj, eegportaddr, trigger_code);  % Turns trigger ON
                wait(eegtriglength);
                io64(ioObj, eegportaddr, 0);     % turns trigger OFF
            end
            if response == 0 % if they haven't made a response yet
                response = 1;
                response_time = GetSecs - start_trial;
            end
        elseif keyCode(rightKey)
            if isEEG
                % Send trigger for trial end
                trigger_code = 400 + trial_info.trials(b,t); % first number: 1 for face, 2 for response, second & third numbers: face type
                io64(ioObj, eegportaddr, trigger_code);  % Turns trigger ON
                wait(eegtriglength);
                io64(ioObj, eegportaddr, 0);     % turns trigger OFF
            end
            if response == 0 % if they haven't made a response yet
                response = 2;
                response_time = GetSecs - start_trial;
            end    
        end
        
    end
    
    if response == 0
        response = NaN;
        response_time = NaN;
    end
    
end


% get response
fprintf(['RT = ' num2str(1000*response_time) 'ms \n']);
accuracy = response == Porientation_order(b,t); % convert button to 1 (hit) or 0 (miss)
if accuracy == 0
    Beeper(250,5,.05);
end
fprintf(['ACC = ' num2str(accuracy) '\n']);
trial_info.practice_acc(b,t) = double(accuracy);
trial_info.practice_rt(b,t) = double(response_time);


