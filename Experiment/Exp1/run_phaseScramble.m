input_folder = 'G:\Jessica\Private\Experiments\2_Experiment2_CFS\EEG_Experiment\Stimuli\Faces';
output_folder = 'G:\Jessica\Private\Experiments\2_Experiment2_CFS\EEG_Experiment\Stimuli\Masks';

face_folders = {'FearfulFemale','FearfulMale','NeutralFemale','NeutralMale'};

iterations = 1:3;

for F = 1:length(face_folders)
    filelist = dir(fullfile(input_folder,face_folders{F},'*.png'));
    for f = 1:length(filelist)
        for i = 1:length(iterations)
            disp(['Phase scrambling ' face_folder{f} ' image ' num2str(f) ', iteration ' num2str(i) '...']);
            phaseScramble(filelist(f).name,fullfile(input_folder,face_folders{F}),output_folder,['scrambled_' num2str(iterations(i)) '_' filelist(f).name]);
        end
    end
end