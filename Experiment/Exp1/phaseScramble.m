function [img_scrambled] = phaseScramble(filename,input_folder,output_folder,output_filename)

% phase scrambles an image (e.g. filename = 'FF01.png', input_folder =
% where the faces are, output_folder = where to save the masks

rand('seed',sum(100*clock));

img = mat2gray(double(imread(fullfile(input_folder,filename))));
random_phase = angle(fft2(rand(size(img,1),size(img,2))));

img_fourier = fft2(img);
Amp = abs(img_fourier);
Phase = angle(img_fourier);
Phase = Phase + random_phase;
img_scrambled = ifft2(Amp .* exp(sqrt(-1) * Phase));

img_scrambled = real(img_scrambled);

imwrite(img_scrambled,fullfile(output_folder,output_filename));

end