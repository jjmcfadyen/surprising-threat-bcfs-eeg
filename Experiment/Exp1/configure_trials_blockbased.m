function [trial_order, image_order, orientation_order, block_order] = configure_trials_blockbased(N,P,M,B,I,O,prob_buffer,block_start)

% N = no. of trials
% P = probability
% M = min_gap
% B = no. of blocks
% I = list of images
% O = orientation order

%%
% FIRST NUMBER                      % SECOND NUMBER
% 1 = Fearful Block, Fearful Face     % 1 = Male
% 2 = Fearful Block, Neutral Face     % 2 = Female
% 3 = Neutral Block, Neutral Face
% 4 = Neutral Block, Fearful Face

if block_start == 1
    block_order = repmat([1 2],1,B/2);
elseif block_start == 2
    block_order = repmat([2 1],1,B/2);
end

trial_order = zeros(B,N);
for b = 1:B
    
    oddball = pseudoRandSpace(P(2)*N,N-prob_buffer,M);
    oddball = [zeros(1,prob_buffer) oddball];
    oddball_idx = find(oddball);
    oddball_randgender = randsample([ones(1,length(oddball_idx/2)) ones(1,length(oddball_idx/2))*2],length(oddball_idx));
    standard_idx = find(~oddball);
    standard_randgender = randsample([ones(1,length(standard_idx/2)) ones(1,length(standard_idx/2))*2],length(standard_idx));
    
    if block_order(b) == 1
        trial_order(b,standard_idx(standard_randgender == 1)) = 31;
        trial_order(b,standard_idx(standard_randgender == 2)) = 32;
        trial_order(b,oddball_idx(oddball_randgender == 1)) = 41;
        trial_order(b,oddball_idx(oddball_randgender == 2)) = 42;
    elseif block_order(b) == 2
        trial_order(b,standard_idx(standard_randgender == 1)) = 11;
        trial_order(b,standard_idx(standard_randgender == 2)) = 12;
        trial_order(b,oddball_idx(oddball_randgender == 1)) = 21;
        trial_order(b,oddball_idx(oddball_randgender == 2)) = 22;
    end
    
end

% codes within 'trial_order'
image_codes =  [21 31;   % Neutral Male
                11 41;   % Fearful Male
                22 32;   % Neutral Female
                12 42];  % Fearful Female
% row/col within 'I'
image_idx =    [1 1;       % Neutral Male
                1 2;       % Fearful Male
                2 1;       % Neutral Female
                2 2];      % Fearful Female

% Randomly select images from entire sample for entire experiment, having
% as many unique face images as possible
image_order = {};
for b = 1:size(trial_order,1)
    for i = 1:size(image_codes,1)
        
        num_images = length(I{image_idx(i,1),image_idx(i,2)});
        num_needed = length(find(trial_order(b,:) == image_codes(i,1) | trial_order(b,:) == image_codes(i,2)));
        
        if num_images >= num_needed
            img_idx = datasample(1:num_images,num_needed,'replace',false);
        end
        
        trial_idx = find(trial_order(b,:) == image_codes(i,1) | trial_order(b,:) == image_codes(i,2));
        
        for t = 1:length(img_idx)

            image_order{b,trial_idx(t)} = I{image_idx(i,1),image_idx(i,2)}{img_idx(t)};

        end
    end
    
    if find(cellfun(@isempty,image_order)) > 0
        error(['Something went wrong in image assignment for block ' num2str(b)]);
    end
    
end

% Which faces are oriented left or right, according to O
all_orientation = [ones(1,O(1)*length(trial_order(:))) ones(1,O(1)*length(trial_order(:)))*2];
orientation_order = randsample(all_orientation,length(all_orientation)); % shuffles
orientation_order = reshape(orientation_order,size(trial_order,1),size(trial_order,2)); % rearranges into B

end