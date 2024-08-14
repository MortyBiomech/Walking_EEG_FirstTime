clc
clear

%% add path to load xdf files
addpath(genpath('C:\Morteza\LSL\xdf-Matlab-master')); % required for loading the XDF files.

file_path = 'C:\Users\morte\OneDrive\Documents\CurrentStudy\sub-P001\ses-S001\eeg';
% file_path = 'C:\Morteza\Analysis\Walking_EEG_Sebastian\datasets\sub-1\ses-S001\eeg';
% file_name = 'sub-1_ses-S001_task-test_reading_GRF_run-001_GRF.xdf';
% GRF = load_xdf(fullfile(file_path, file_name));

% file_name = 'sub-1_ses-S001_task-test_12August_run-001_eeg.xdf';
file_name = 'sub-P001_ses-S001_task-Default_run-001_eeg.xdf';

Data = load_xdf(fullfile(file_path, file_name));



%% Plot 8 channels separately
tiledlayout(4,2)
nexttile
plot(GRF{1,1}.time_stamps, GRF{1,1}.time_series(1, :))

nexttile
plot(GRF{1,1}.time_stamps, GRF{1,1}.time_series(2, :), 'r')

nexttile
plot(GRF{1,1}.time_stamps, GRF{1,1}.time_series(3, :), 'r')

nexttile
plot(GRF{1,1}.time_stamps, GRF{1,1}.time_series(4, :))

nexttile
plot(GRF{1,1}.time_stamps, GRF{1,1}.time_series(5, :))

nexttile
plot(GRF{1,1}.time_stamps, GRF{1,1}.time_series(6, :), 'r')

nexttile
plot(GRF{1,1}.time_stamps, GRF{1,1}.time_series(7, :), 'r')

nexttile
plot(GRF{1,1}.time_stamps, GRF{1,1}.time_series(8, :))

%% plot results (Left and Right GRF)
Left_leg_indx = [2 3 6 7];
Riht_leg_indx = [1 4 5 8];

figure()
GRF = Data{1, 3};
plot(GRF.time_stamps, sum(GRF.time_series(Left_leg_indx, :), 1), 'r')
hold on
plot(GRF.time_stamps, sum(GRF.time_series(Riht_leg_indx, :), 1), 'b')


%% plot GRF and EMG 
tiledlayout(3, 1)
nexttile
plot(GRF_and_EMG{1,1}.time_stamps, sum(GRF_and_EMG{1,1}.time_series(Left_leg_indx, :), 1), 'r')
hold on
plot(GRF_and_EMG{1,1}.time_stamps, sum(GRF_and_EMG{1,1}.time_series(Riht_leg_indx, :), 1), 'b')

nexttile
plot(GRF_and_EMG{1,2}.time_stamps, GRF_and_EMG{1,2}.time_series(2,:), 'b')

nexttile
plot(GRF_and_EMG{1,2}.time_stamps, GRF_and_EMG{1,2}.time_series(3,:), 'r')
