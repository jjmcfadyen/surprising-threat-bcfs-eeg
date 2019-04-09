function [tex, mon_load] = load_stimuli(text_pos,trial_order,image_order,mondrians,window)

DrawFormattedText(window,['Loading images... ' num2str(0) '%'],text_pos,'center');
Screen('Flip',window);

% Faces
total_trials = length(trial_order(:));
counter = 0;
tex = [];
for b = 1:size(image_order,1)
    for t = 1:size(image_order,2)
        counter = counter + 1;
        progress = round((counter/total_trials)*100);
        DrawFormattedText(window,['Loading images... ' num2str(progress) '%'],text_pos,'center');
        Screen('Flip',window);
        tex(b,t) = Screen('MakeTexture',window,imread(char(image_order{b,t})));
    end
end

% Masks
mon_load = [];
for m = 1:length(mondrians)
    mon_load(m) = Screen('MakeTexture',window,imread(['Stimuli/Masks/' mondrians{m}]));
    DrawFormattedText(window,['Loading masks... ' num2str(floor(m/length(mondrians)*100)) '%'],text_pos,'center');
    Screen('Flip',window);
end
