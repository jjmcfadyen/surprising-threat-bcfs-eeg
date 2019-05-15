%% List of bad/noisy channels from Experiment 2
% includes study notes

bad_channels = {};
bad_channels{1} = {'Fpz'};
bad_channels{2} = {'Fpz'};
bad_channels{3} = {'Fpz','Fz'};
bad_channels{4} = {'Fpz','Fz'};
bad_channels{5} = {'Fpz','Fz'};
bad_channels{6} = {'Fpz','Fz'};
bad_channels{7} = {'Fpz','Fz'};
bad_channels{8} = {'Fpz','Fz'};
bad_channels{9} = {'Fpz','Fz'};
bad_channels{10} = {'Fpz','Fz','AFz','P2'};
bad_channels{11} = {'Fpz','Fz','AFz'};
bad_channels{12} = {'Fpz','Fz','AFz'};
bad_channels{13} = {'Fpz','Fz','AFz'};
bad_channels{14} = {'Fpz','Fz','AFz','Pz'};
bad_channels{15} = {'Fpz','Fz','AFz'};
bad_channels{16} = {'Fpz','Fz','AFz'};
bad_channels{17} = {'Fpz','Fz','AFz','P2'};
bad_channels{18} = {'Fpz','Fz','AFz'};
bad_channels{19} = {'Fpz','Fz','AFz'};
bad_channels{20} = {'Fpz','Fz','AFz'};
bad_channels{21} = {'Fpz','Fz','AFz'};
bad_channels{22} = {'Fpz','Fz','AFz'};
bad_channels{23} = {'Fpz','Fz','AFz'};
bad_channels{24} = {'Fpz','Fz','AFz'};
bad_channels{25} = {'Fpz','Fz','AFz'};
bad_channels{26} = {'Fpz','Fz','AFz'}; % F5 had weird pulse... P2 had heart rate artefact. Strong drift
bad_channels{27} = {'Fpz','Fz'}; % EXG1 = 'AFz', EXG2 = 'F2'
bad_channels{28} = {'Fpz','Fz'}; % EXG1 = 'AFz', EXG2 = 'F2'
bad_channels{29} = {'Fpz','Fz','AFz','F2'};
bad_channels{30} = {'Fpz','Fz'}; % EXG1 = 'AFz', EXG2 = 'F2'
bad_channels{31} = {'Fpz','Fz','C2','FCz','P10'}; % EXG1 = 'AFz', EXG2 = 'F2'. T8 noisy towards end
bad_channels{32} = {'Fpz','Fz','P10'}; % EXG1 = 'AFz', EXG2 = 'F2'
bad_channels{33} = {'Fpz','Fz','P10'}; % EXG1 = 'AFz', EXG2 = 'F2' (MIXED UP EXG2 and EXG8!)

swapped_channels = cell(1,33);
swapped_channels{27} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{28} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{30} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{31} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{32} = {'Nz','AFz'; 'SNz','F2'};
swapped_channels{33} = {'Nz','AFz'; 'M2','F2'};
