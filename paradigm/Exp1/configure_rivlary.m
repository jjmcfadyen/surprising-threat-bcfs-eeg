%% Configure Rivalry

% Draw centre line
centre_line = CenterRectOnPointd([0 0 3 screenYpixels],xCentre,yCentre);
Screen('FillRect', window, [0 0 0], centre_line);
Screen('Flip',window,[],1);

vertical_pos = yCentre;

% Draw squares and move with left/right button presses
end_conf = 0;
while end_conf == 0
    
    % Draw centre line
    centre_line = CenterRectOnPointd([0 0 3 screenYpixels],xCentre,yCentre);
    Screen('FillRect', window, [0 0 0], centre_line);
    
    leftRect_outer = CenterRectOnPointd([0 0 stim_size*border_width stim_size*border_width],xCentre-from_centre, vertical_pos);
    rightRect_outer = CenterRectOnPointd([0 0 stim_size*border_width stim_size*border_width],xCentre+from_centre, vertical_pos);
    leftRect_inner = CenterRectOnPointd([0 0 stim_size*1.25 stim_size*1.25],xCentre-from_centre, vertical_pos);
    rightRect_inner = CenterRectOnPointd([0 0 stim_size*1.25 stim_size*1.25],xCentre+from_centre, vertical_pos);
    % draw squares
    Screen('FillRect', window, [0 0 0], leftRect_outer);  % black border left
    Screen('FillRect', window, [0 0 0], rightRect_outer); % black border right
    Screen('FillRect', window, [1 1 1], leftRect_inner);  % red square left
    Screen('FillRect', window, [1 1 1], rightRect_inner); % blue square right
    % draw fixation crosses
    DrawFormattedText(window,'+','center','center',0,[],[],[],[],[],leftRect_inner);
    DrawFormattedText(window,'+','center','center',0,[],[],[],[],[],rightRect_inner);
    Screen('Flip',window);
    
    [keyIsDown,secs,keyCode] = PsychHID('KbCheck',[],[]);
    if find(keyCode) == leftKey
        from_centre = from_centre + 1;
    elseif find(keyCode) == rightKey
        from_centre = from_centre - 1;
    elseif find(keyCode) == upKey
        vertical_pos = vertical_pos - 1;
    elseif find(keyCode) == downKey
        vertical_pos = vertical_pos + 1;    
    elseif find(keyCode) == spaceKey
        break
    end
end

% Update 'from_centre' variable
left_centre = xCentre-from_centre;
right_centre = xCentre+from_centre;
