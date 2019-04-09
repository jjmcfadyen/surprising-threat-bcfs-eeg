%% Run Trials

% Prepare stimuli indices
if ~face_flicker
        inc = 1/(stim_dur*monRate);
else
    T = 1:((stim_dur+this_ITI)*monRate); % no. of frames per trial
    t = 1:length(T)/(face_flicker*stim_dur); % no. of frames per fade in & out
    w = (sin(2*pi*(1/length(t))*t+(1.5*pi))+1)/2;
    

    if abs(min_face_contrast - max_face_contrast) > 0
        face_fade = min_face_contrast:(max_face_contrast-min_face_contrast)/(face_flicker*stim_dur):max_face_contrast;
        face_fade = face_fade(2:end);
    else face_fade = ones(1,ceil(face_flicker*stim_dur))*min_face_contrast;
    end
    face_idx = [];
    for i = 1:face_flicker*stim_dur
        face_idx = [face_idx w*face_fade(i)];
    end
    if this_ITI > 0
        face_idx = [face_idx zeros(1,ceil((monRate/this_ITI)+1))]; % add extra to end
    else face_idx = [face_idx 1];
    end
end

% Set up presentation boxes
if dominant_eye == 1
    this_dom = left_centre;
    this_sub = right_centre;
elseif dominant_eye == 2
    this_dom = right_centre;
    this_sub = left_centre;
end
domRect_inner = CenterRectOnPointd([0 0 stim_size*1.25 stim_size*1.25],this_dom, vertical_pos);
subRect_inner_face = CenterRectOnPointd([0 0 stim_size stim_size],this_sub, vertical_pos);
domRect_outer = CenterRectOnPointd([0 0 stim_size*border_width stim_size*border_width],this_dom, vertical_pos);
subRect_outer = CenterRectOnPointd([0 0 stim_size*border_width stim_size*border_width],this_sub, vertical_pos);
domRect_inner_face = CenterRectOnPointd([0 0 stim_size stim_size],this_dom, vertical_pos);

leftRect = CenterRectOnPointd([0 0 screenXpixels/2 screenYpixels],left_centre, yCentre);
rightRect = CenterRectOnPointd([0 0 screenXpixels/2 screenYpixels],right_centre, yCentre);

%% Present stimuli

if practice_now
    these_trials = size(trial_info.practice_trials,2);
else these_trials = size(trial_info.trials,2);
end

mask_idx = 0:monRate/mask_flicker:(monRate*(stim_dur+max(ITI)))+100;

refresh_log = [];
for t = 1:these_trials
    
    fprintf(['TRIAL ' num2str(t) ' of ' num2str(size(trial_order,2)) '\n']);
    this_ITI = randsample(ITI,1);
    
    % Check if control condition
    if ~practice_now
        if mask_control
            this_mask_control = trial_info.mask_control(b,t);
        else
            this_mask_control = 0;
        end
    else this_mask_control = 0;
    end

    % Which orientation?
    if ~practice_now
        if orientation_order(b,t) == 1
            rotation =  -rotation_degree; % left
        elseif orientation_order(b,t) == 2
            rotation = rotation_degree; % right
        end
    elseif practice_now
        if Porientation_order(b,t) == 1
            rotation =  -rotation_degree; % left
        elseif Porientation_order(b,t) == 2
            rotation = rotation_degree; % right
        end
    end
    
    % Which trigger?
    if isEEG
        if ~practice_now
            trigger_code = trial_info.eegtrig(b,t);
        elseif practice_now
            trigger_code = trial_info.practice_eegtrig(b,t);
        end
    end

    %% Get presentation structure for trial
    
    mondrian_counter = 1;
    face_counter = 1; % when to change contrast
    curr_mask = 1;
    response = 0;

    if isEEG
        start_trigger = 0;
    end
    
    frame_counter = 0;
    start_trial = GetSecs; % start reaction time timer
    FC = linspace(min_face_contrast,max_face_contrast,(monRate*stim_dur)+1);
    while frame_counter < monRate*stim_dur
        
        frame_counter = frame_counter + 1;
        
        % draw boxes
        Screen('FillRect', window, [0 0 0], domRect_outer);
        Screen('FillRect', window, [0 0 0], subRect_outer);
        
        % select mondrian
        if frame_counter > mask_idx(mondrian_counter)
            mondrian_counter = mondrian_counter + 1;
            mask_options = find(1:length(masks)~=curr_mask);
            curr_mask = mask_options(randi(length(mask_options)));
            curr_mask_rotation = datasample(0:90:360,1);
        end

        % draw face
        if ~face_flicker

%             face_contrast = min(max_face_contrast,min_face_contrast+(inc*frame_counter));
            face_contrast = FC(frame_counter)
            Screen('DrawTexture',window,tex(b,t),[],subRect_inner_face,rotation,[],face_contrast);
            if this_mask_control
                Screen('DrawTexture',window,tex(b,t),[],domRect_inner_face,rotation,[],face_contrast);
            end
            
        elseif face_flicker
            
            % draw face
            if dummy_tex
                Screen('FillRect', window, [1 1 1]*face_idx(frame_counter), subRect_inner_face);
                if this_mask_control
                    Screen('FillRect', window, [1 1 1]*face_idx(frame_counter), domRect_inner_face);
                end
            else
                Screen('DrawTexture',window,tex(b,t),[],subRect_inner_face,rotation,[],face_idx(frame_counter));
                if this_mask_control
                    Screen('DrawTexture',window,tex(b,t),[],domRect_inner_face,rotation,[],face_idx(frame_counter));
                end
            end
            
        end

        % draw mondrian
        if dummy_tex
            if mod(mondrian_counter,2)
                Screen('FillRect', window, [0 0 0], domRect_outer);
            else Screen('FillRect', window, [1 1 1], domRect_outer);
            end
        else
            if ~this_mask_control
                Screen('DrawTexture',window,mon_load( curr_mask),[],domRect_inner,curr_mask_rotation,[],mondrian_contrast);
            end
        end

        % draw fixation crosses
        DrawFormattedText(window,'+',left_centre-(fontsize/2),vertical_pos-(fontsize/2),0,[],[],[],[],[],leftRect);
        DrawFormattedText(window,'+',right_centre-(fontsize/2),vertical_pos-(fontsize/2),0,[],[],[],[],[],rightRect);

        % flip
        [VBLTimestamp StimulusOnsetTime FlipTimestamp] = Screen('Flip',window);
        
        if isEEG
            if ~start_trigger
                % Send trigger for trial start
                if isEEG
                    io64(ioObj, eegportaddr, trigger_code);  % Turns trigger ON
                    WaitSecs(eegtriglength);
                    io64(ioObj, eegportaddr, 0);     % turns trigger OFF
                end
                start_trigger = 1;
            end
        end

        % check for response
        [keyIsDown,secs,keyCode] = KbCheck;
        if keyCode(escapeKey)
            sca;
            return;
        elseif keyCode(leftKey)
            if response == 0 % if they haven't made a response yet
                % Send trigger for response
                if isEEG
                    io64(ioObj, eegportaddr, trigger_code+20);  % Turns trigger ON
                    WaitSecs(eegtriglength);
                    io64(ioObj, eegportaddr, 0);     % turns trigger OFF
                end
                response = 1;
                response_time = GetSecs - start_trial;
            end
        elseif keyCode(rightKey)
            if response == 0 % if they haven't made a response yet
                 % Send trigger for response
                if isEEG
                    io64(ioObj, eegportaddr, trigger_code+20);  % Turns trigger ON
                    WaitSecs(eegtriglength);
                    io64(ioObj, eegportaddr, 0);     % turns trigger OFF
                end
                response = 2;
                response_time = GetSecs - start_trial;
            end
        end
        
        if face_flicker
            refresh_log = [refresh_log; t frame_counter VBLTimestamp StimulusOnsetTime FlipTimestamp face_idx(frame_counter)];
        else refresh_log = [refresh_log; t frame_counter VBLTimestamp StimulusOnsetTime FlipTimestamp face_contrast];
        end
        
    end
    
    %% Masked ITI
    if this_ITI > 0
        
        % Get presentation structure for trial
        mondrian_counter = 1;
        curr_mask = 1;
        
        % Send trigger for face offset
        if isEEG
            io64(ioObj, eegportaddr, 99);  % Turns trigger ON
            WaitSecs(eegtriglength);
            io64(ioObj, eegportaddr, 0);     % turns trigger OFF
        end
            
        frame_counter = 0;
        while frame_counter <= monRate*this_ITI

            frame_counter = frame_counter+1;
            
            % draw boxes
            Screen('FillRect', window, [0 0 0], domRect_outer);
            Screen('FillRect', window, [0 0 0], subRect_outer);
            
            if ~this_mask_control
                % select mondrian
                if frame_counter > mask_idx(mondrian_counter)
                    mondrian_counter = mondrian_counter + 1;
                    mask_options = find(1:length(masks)~=curr_mask);
                    curr_mask = mask_options(randi(length(mask_options)));
                    curr_mask_rotation = datasample(0:90:360,1);
                end

                % draw mondrian
                if dummy_tex
                    if mod(mondrian_counter,2)
                        Screen('FillRect', window, [0 0 0], domRect_outer);
                    else Screen('FillRect', window, [1 1 1], domRect_outer);
                    end
                else
                    Screen('DrawTexture',window,mon_load( curr_mask),[],domRect_inner,curr_mask_rotation,[],mondrian_contrast);
                end

                % draw fixation crosses
                DrawFormattedText(window,'+',left_centre-(fontsize/2),vertical_pos-(fontsize/2),0,[],[],[],[],[],leftRect);
                DrawFormattedText(window,'+',right_centre-(fontsize/2),vertical_pos-(fontsize/2),0,[],[],[],[],[],rightRect);
            end

            % flip
            Screen('Flip',window);

        end
    end

    
    if response == 0
        response = NaN;
        response_time = NaN;
    end

    % get response
    fprintf(['RT = ' num2str(1000*response_time) 'ms \n']);
    if ~practice_now
        accuracy = response == orientation_order(b,t); % convert button to 1 (hit) or 0 (miss)
        fprintf(['ACC = ' num2str(accuracy) '\n']);
        trial_info.acc(b,t) = double(accuracy);
        trial_info.rt(b,t) = double(response_time);
    elseif practice_now
        accuracy = response == Porientation_order(b,t); % convert button to 1 (hit) or 0 (miss)
        fprintf(['ACC = ' num2str(accuracy) '\n']);
        trial_info.practice_acc(b,t) = double(accuracy);
        trial_info.practice_rt(b,t) = double(response_time); 
    end

end
