clc
clear

%% Change these paths with respect to your system
main_data_path = 'C:\Morteza\MyProjects\Walking_EEG_FirstTime\data';

study_path1 = 'C:\Morteza\MyProjects\Walking_EEG_FirstTime\data\Walking_Different_Speeds';
study_path2 = 'C:\Morteza\MyProjects\Walking_EEG_FirstTime\data\Walking_Different_Weights';
study_path = study_path1;

if strcmp(study_path, study_path1)
    study_identifier = 1;
else
    study_identifier = 2;
end

processing_path = 'C:\Morteza\MyProjects\Walking_EEG_FirstTime\Code\Matlab\data_processing\';
addpath(genpath('C:\Morteza\MyProjects\Walking_EEG_FirstTime'))

%% Add important paths
addpath('C:\Morteza\Toolboxes\Fieldtrip\fieldtrip-20231127')
addpath('C:\Morteza\Toolboxes\Fieldtrip\fieldtrip-20231127\fileio')

%% All signals from all sessions concatenated (it takes time!)
subject_id = 2;
rawdata_path = [main_data_path, filesep, '0_source_data\'];
output = extract_streams(subject_id, rawdata_path, study_identifier);

All_EEG_time = output.All_EEG_time;
All_GRF = output.All_GRF;
All_GRF_time = output.All_GRF_time;

Left_leg_indx = [2 3 6 7];
Right_leg_indx = [1 4 5 8];

All_GRF_Left  = sum(All_GRF(Left_leg_indx, :), 1);
All_GRF_Right = sum(All_GRF(Right_leg_indx, :), 1);


%% Initialize EEGLAB 
if ~exist('ALLCOM','var')
	eeglab;
end

% initialize fieldtrip without adding alternative files to path 
% assuming FT is on your path already or is added via EEGlab plugin manager
global ft_default
ft_default.toolbox.signal = 'matlab';  % can be 'compat' or 'matlab'
ft_default.toolbox.stats  = 'matlab';
ft_default.toolbox.image  = 'matlab';
ft_defaults % this sets up the FieldTrip path


%% [OPTIONAL] check the .xdf data to explore the structure
ftPath = fileparts(which('ft_defaults')); 
addpath(fullfile(ftPath, 'external','xdf')); 
xdfPath = [rawdata_path, 'sub-', num2str(subject_id), filesep, ...
    'ses-S001\eeg\sub-', num2str(subject_id), ...
    '_ses-S001_task-WalkingDifferentSpeeds_run-001_eeg.xdf']; % enter full path to your .xdf file 

% load .xdf data to check what is in there
streams         = load_xdf(xdfPath);
streamnames     = cellfun(@(x) x.info.name, streams, 'UniformOutput', 0)' % will display names of streams contained in .xdf

% display names of all channels in the .xdf data
for Si = 1:numel(streamnames)
    if isfield( streams{Si}.info.desc, 'channels')
        channelnames    = cellfun(@(x) x.label, streams{Si}.info.desc.channels.channel, 'UniformOutput', 0)'
    end
end


%% [OPTIONAL] enter metadata about the data set, data modalities, and participants

% information about the eeg recording system 
% will be saved in BIDS-folder/sub-XX/[ses-XX]/eeg/*_eeg.json and *_coordsystem.json
% see "https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#:~:text=MAY%20be%20specified.-,Sidecar%20JSON%20(*_eeg.json),-Generic%20fields%20MUST"
% and "https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#:~:text=after%20the%20recording.-,Coordinate%20System%20JSON,-(*_coordsystem.json"
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
eegInfo     = [];
% eegInfo.coordsystem.EEGCoordinateSystem     = 'enter the name of your coordinate system'; % only needed when you share eloc
% eegInfo.coordsystem.EEGCoordinateUnits      = 'enter the unit of your coordinate system'; % only needed when you share eloc
% eegInfo.coordsystem.EEGCoordinateSystemDescription = 'enter description of your coordinate system'; % only needed when you share eloc
eegInfo.eeg.SamplingFrequency               = 500; % nominal sampling frequency  
                               


%% Import data

% iterate over sessions to import xdf files one by one
% full path to your study folder  
studyFolder   = study_path(1:end-1);     

% replace with names of your sessions (if there are no multiple sessions, 
% remove confg.ses and the session loop in the following)
sessionNames  = {'S001', 'S002', 'S003', 'S004'};             


% loop over sessions 
for session = 1:length(sessionNames)
    
    % reset for each loop
    config = [];  
    % required, replace with the folder where you want to store your bids data
    config.bids_target_folder = [study_path, '1_BIDS_data'];  
    
    % required, replace with your xdf file full path
    config.filename  = ...
        fullfile([rawdata_path,'sub-', num2str(subject_id),'\ses-', ...
        sessionNames{session}, '\eeg\sub-', num2str(subject_id), ...
        '_ses-', sessionNames{session}, '_task-Default_run-001_eeg.xdf']);    
    
    % optional, if you have electrode location file, replace with the full 
    % path to the file.
    % config.eeg.chanloc = fullfile([processing_path ,'chanlocs.ced']);
    
    % optional, replace with your task name
    config.task                   = 'KneeSwingingwithExo';  
    config.subject                = subject_id;  % required
    config.session                = sessionNames{session}; % optional
    config.overwrite              = 'on'; % optional
    
    % required, replace with the unique keyword in your eeg stream in 
    % the .xdf file
    config.eeg.stream_name        = 'LiveAmpSN-102108-1125'; 
    
    %------------------------------------------------------------------
    
    % config.motion.streams{1}.xdfname            = 'YourStreamNameInXDF'; % replace with name of the stream corresponding to the first tracking system
    % config.motion.streams{1}.bidsname           = tracking_systems{1};  % a comprehensible name to represent the tracking system 
    % config.motion.streams{1}.tracked_points     = 'headRigid';          % name of the point that is being tracked in the tracking system, the keyword has to be containted in the channel name (see "bemobil_bids_motionconvert")
    % config.motion.streams{1}.tracked_points_anat= 'head';               % example of how the tracked point can be renamed to body part name for metadata
    % 
    % % names of position and quaternion channels in each stream
    % config.motion.streams{1}.positions.channel_names    = {'headRigid_Rigid_headRigid_X';  'headRigid_Rigid_headRigid_Y' ; 'headRigid_Rigid_headRigid_Z' };
    % config.motion.streams{1}.quaternions.channel_names  = {'headRigid_Rigid_headRigid_quat_W';'headRigid_Rigid_headRigid_quat_Z';...
    %                                                        'headRigid_Rigid_headRigid_quat_X';'headRigid_Rigid_headRigid_quat_Y'};
    % 
    % config.motion.streams{2}.xdfname            = 'YourStreamNameInXDF2';
    % config.motion.streams{2}.bidsname           = tracking_systems{2};
    % config.motion.streams{2}.tracked_points     = {'Rigid1', 'Rigid2', 'Rigid3', 'Rigid4'}; % example when there are multiple points tracked by the system
    % config.motion.streams{2}.positions.channel_names = {'Rigid1_X', 'Rigid2_X', 'Rigid3_X', 'Rigid4_X';... % each column is one tracked point and rows are different coordinates
    %                                                     'Rigid1_Y', 'Rigid2_Y', 'Rigid3_Y', 'Rigid4_Y';...
    %                                                     'Rigid1_Z', 'Rigid2_Z', 'Rigid3_Z', 'Rigid4_Z'};
    % config.motion.streams{2}.quaternions.channel_names = {'Rigid1_A', 'Rigid2_A', 'Rigid3_A', 'Rigid4_A';...
    %                                                     'Rigid1_B', 'Rigid2_B', 'Rigid3_B', 'Rigid4_B';...
    %                                                     'Rigid1_C', 'Rigid2_C', 'Rigid3_C', 'Rigid4_C'; ...
    %                                                     'Rigid1_D', 'Rigid2_D', 'Rigid3_D', 'Rigid4_D'};
    % 

    % bemobil_xdf2bids(config, ...
    %     'general_metadata', generalInfo,...
    %     'participant_metadata', subjectInfo,...
    %     'motion_metadata', motionInfo, ...
    %     'eeg_metadata', eegInfo);

    bemobil_xdf2bids(config, ...
        'eeg_metadata', eegInfo);
end


fclose all;


%% configuration for bemobil bids2set
%----------------------------------------------------------------------
config.set_folder               = fullfile(studyFolder,'2_raw-EEGLAB');
config.session_names            = sessionNames;

% config.other_data_types = {'motion'};  % specify which other data type than eeg is there (only 'motion' and 'physio' supported atm)
bemobil_bids2set(config);


%% Configuration of the BeMoBIL pipeline
bemobil_config = BeMoBIL_Configuration(study_path);

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = 5; 

% set to 1 if all files should be computed, independently of whether they are present on disk or not
force_recompute = 0; 



%% processing loop

for subject = subjects
    
    %% prepare filepaths and check if already done
    
	disp(['Subject #' num2str(subject)]);
    
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = [];  CURRENTSET=[]; EEG=[]; EEG_interp_avref = []; EEG_single_subject_final = [];
	
	input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
	output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_analysis_folder bemobil_config.filename_prefix num2str(subject)];
	
	try
		% load completely processed file
		EEG_single_subject_final = ...
            pop_loadset('filename', [ bemobil_config.filename_prefix num2str(subject)...
			'_' bemobil_config.single_subject_cleaned_ICA_filename], 'filepath', output_filepath);
    catch
        disp('...failed. Computing now.')
    end

	
	if ~force_recompute && exist('EEG_single_subject_final','var') ...
            && ~isempty(EEG_single_subject_final)
		clear EEG_single_subject_final
		disp('Subject is completely preprocessed already.')
		continue
    end

	%% load data in EEGLAB .set structure
    % make sure the data is stored in double precision, large datafiles are
    % supported, no memory mapped objects are used but data is processed 
    % locally, and two files are used for storing sets (.set and .fdt)
	try 
        pop_editoptions('option_saveversion6', 0, 'option_single', 0, ...
            'option_memmapdata', 0, 'option_savetwofiles', 1, ...
            'option_storedisk', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!'); 
    end
    
    % load files that were created from xdf to BIDS to EEGLAB
    EEG = pop_loadset('filename', ...
        [ bemobil_config.filename_prefix num2str(subject) '_' ...
        bemobil_config.merged_filename],'filepath',input_filepath);
    

    %% reselect EEG channels and remove ACC channels
    EEG = pop_select(EEG, 'rmchannel', [65, 66, 67]);

    %% Load channels location file
    EEG = pop_chanedit(EEG,'load', [processing_path, 'chanlocs.ced']); 

    %% Define and Add events
    % Compute latency values
    
    % start_beep event
    start_beep = find(diff(All_Experiment(6, :)) == 1);
    start_beep(1:6) = [];
    start_beep_time_Expdata = All_Experiment_time(start_beep);
    start_beep_indx_EEG = ...
        knnsearch(All_EEG_time', start_beep_time_Expdata');
    start_beep_latency_EEG = EEG.times(start_beep_indx_EEG); % time unit: milisecond
    
    % pressure_change event
    pressure_change_time_Expdata = All_Experiment_time(start_beep) - 2;
    pressure_change_indx_EEG = ...
        knnsearch(All_EEG_time', pressure_change_time_Expdata');
    pressure_change_latency_EEG = EEG.times(pressure_change_indx_EEG); % time unit: milisecond
    
    % finish_beep event
    finish_beep = find(diff(All_Experiment(6, :)) == -1);
    finish_beep(1:6) = [];
    finish_beep_time_Expdata = All_Experiment_time(finish_beep);
    finish_beep_indx_EEG = ...
        knnsearch(All_EEG_time', finish_beep_time_Expdata');
    finish_beep_latency_EEG = EEG.times(finish_beep_indx_EEG); % time unit: milisecond
    
    % Trial_Start event (500ms before pressure_change for baseline removal)
    Trial_Start_time_Expdata = All_Experiment_time(start_beep) - 2.5;
    Trial_Start_indx_EEG = ...
        knnsearch(All_EEG_time', Trial_Start_time_Expdata');
    Trial_Start_latency_EEG = EEG.times(Trial_Start_indx_EEG); % time unit: milisecond

    % Trial_End event (500ms after finish_beep as subject were asked to
    % stop moving their knee and evaluate the task after hearing the
    % finish-beep sound (double-beeps).
    Trial_End_time_Expdata = All_Experiment_time(finish_beep) + 0.5;
    Trial_End_indx_EEG = ...
        knnsearch(All_EEG_time', Trial_End_time_Expdata');
    Trial_End_latency_EEG = EEG.times(Trial_End_indx_EEG); % time unit: milisecond



    % Preallocate the memory for score_change_latency_EEG
    score_press_latency_EEG = zeros(1, numel(start_beep));
    
    % Create an interpolation function for EEG times
    eeg_time_interpolant = griddedInterpolant(All_EEG_time, 1:numel(All_EEG_time), 'nearest');
    
    for i = 1:numel(start_beep)-1
        % Get the segment of Expdata_temp and Exptime_temp
        Expdata_temp = All_Experiment(4, finish_beep(i):start_beep(i+1));
        Exptime_temp = All_Experiment_time(1, finish_beep(i):start_beep(i+1));
    
        % Find the index where score change occurs
        score_press_indx = find(abs(Expdata_temp - Expdata_temp(1)) >= 1, 1);
    
        % This should be changed (Add the moment of pressing score buttons)
        if ~isempty(score_press_indx) % score changed
            score_press_time_Expdata = Exptime_temp(score_press_indx);
        else % score didn't change, assume 2 seconds after finish_beep
            score_press_time_Expdata = All_Experiment_time(finish_beep(i)) + 2;
        end
    
        % Use the interpolation function to find the closest EEG index
        score_press_indx_EEG = round(eeg_time_interpolant(score_press_time_Expdata));
    
        % Calculate latency
        score_press_latency_EEG(i) = EEG.times(score_press_indx_EEG);
    end
    
    % for the last trial
    % Get the segment of Expdata_temp and Exptime_temp
    Expdata_temp = All_Experiment(4, finish_beep(end):end);
    Exptime_temp = All_Experiment_time(1, finish_beep(end):end);
    
    % Find the index where score change occurs
    score_press_indx = find(abs(Expdata_temp - Expdata_temp(1)) >= 1, 1);
    
    if ~isempty(score_press_indx) % score was changed
        score_press_time_Expdata = Exptime_temp(score_press_indx);
    else % score wasn't change, assume 2 seconds after finish_beep
        score_press_time_Expdata = All_Experiment_time(finish_beep(end)) + 2;
    end
    
    % Use the interpolation function to find the closest EEG index
    score_press_indx_EEG = round(eeg_time_interpolant(score_press_time_Expdata));
    
    % Calculate latency
    score_press_latency_EEG(end) = EEG.times(score_press_indx_EEG);
    
    
    %% Define and Add Events
    % Import Trials information
    Trials = cell(1, numel(start_beep));
    for i = 1:numel(Trials)
        Trials{1, i}.Pressure = All_Experiment(3, start_beep(1, i));
        if i ~= numel(Trials)
            Trials{1, i}.Score = All_Experiment(4, start_beep(1, i+1));
        else
            Trials{1, i}.Score = All_Experiment(4, end);
        end
    end
    
    % Add pressure change events
    type = repmat({'PC_Pressure_Change'}, 1, numel(pressure_change_latency_EEG));
    latency = pressure_change_latency_EEG;
    desc1 = cell(1, numel(pressure_change_latency_EEG)); % {P(i-1), P(i), Trial}
    for i = 1:length(desc1)
        if i~=1
            desc1{1, i} = {Trials{1, i-1}.Pressure, Trials{1, i}.Pressure, i};
        else
            indx_temp = knnsearch(All_Experiment_time', pressure_change_time_Expdata(1) - 0.1);
            desc1{1, 1} = {All_Experiment(3, indx_temp), Trials{1, i}.Pressure, i};
        end
    end
    desc = desc1;
    
    % Add start beep events
    type = cat(2, type, repmat({'SB_Start_Beep'}, 1, numel(start_beep)));
    latency = cat(2, latency, start_beep_latency_EEG);
    desc = cat(2, desc, desc1);
    
    % Add finish beep events
    type = cat(2, type, repmat({'FB_Finish_Beep'}, 1, numel(finish_beep)));
    latency = cat(2, latency, finish_beep_latency_EEG);
    desc = cat(2, desc, desc1);

    % Add score change events
    type = cat(2, type, repmat({'SP_Score_Press'}, 1, numel(score_press_latency_EEG)));
    latency = cat(2, latency, score_press_latency_EEG);
    desc2 = cell(1, numel(score_press_latency_EEG)); % {P(i-1), P(i), Trial}
    for i = 1:length(desc2)
        if i~=1
            desc2{1, i} = {Trials{1, i-1}.Pressure, Trials{1, i}.Pressure,Trials{1, i-1}.Score,Trials{1, i}.Score, i};
        else
            indx_temp = knnsearch(All_Experiment_time', pressure_change_time_Expdata(1) - 0.1);
            desc2{1, 1} = {All_Experiment(3, indx_temp), Trials{1, i}.Pressure, All_Experiment(4, indx_temp), Trials{1, i}.Score, i};
        end
    end
    desc = cat(2, desc, desc2);

    % Add Trial Start events
    type = cat(2, type, repmat({'TS_Trial_Start'}, 1, numel(Trial_Start_latency_EEG)));
    latency = cat(2, latency, Trial_Start_latency_EEG);
    desc = cat(2, desc, desc1);

    % Add Trial End events
    type = cat(2, type, repmat({'TE_Trial_End'}, 1, numel(Trial_End_latency_EEG)));
    latency = cat(2, latency, Trial_End_latency_EEG);
    desc = cat(2, desc, desc1);

    

    %% Write the TS_PC_SB_FB_TE_SP_event.txt file
    % TS: Trial Start
    % PC: Pressure Change
    % SB: Sstart Beep
    % FB: Finish Beep
    % TE: Trial End
    % SP: Score Press

    folder = [processing_path, 'Events', ...
        filesep, 'sub-', num2str(subject_id)];

    % Ensure the folder exists, if not, create it
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    % File name
    filename = fullfile(folder, 'events.txt');
    
    % Open the file for writing
    fileID = fopen(filename, 'w');
    
    % Check if the file was opened successfully
    if fileID == -1
        error('Cannot open file for writing: %s', filename);
    end
    
    % Write the header
    fprintf(fileID, 'type\tlatency\tdesc\n');
    
    % Write the data
    for i = 1:numel(type)
        % Convert the nested cell array in desc to a string with underline separator
        desc_str = strjoin(cellfun(@num2str, desc{i}, 'UniformOutput', false), '_');
        fprintf(fileID, '%s\t%d\t%s\n', type{i}, latency(i), desc_str);
    end
    
    % Close the file
    fclose(fileID);
    
    % Notify the user
    fprintf('File saved successfully: \n%s\n', filename);
    
    
    %% Add Events to the EEG file 
    [EEG, eventnumbers] = pop_importevent(EEG, 'event', ...
              filename, 'fields', {'type', 'latency','desc' }, ...
              'append', 'no', 'align', NaN, 'skipline', 1, 'timeunit', 1E-3);


    %% individual EEG processing to remove non-exp segments
    % it is stongly recommended to remove these segments because they may contain strong artifacts that confuse channel
    % detection and AMICA
    
    % this example removes everything before the first and after the last event with a buffer of 1 second
    
    allevents = {EEG.event.type}';
    
    removeindices = zeros(numel(start_beep)+1 ,2);
    % remove from start to first event
    removeindices(1, :) = [0 EEG.event(1).latency-50]; 

    % add more removeIndices here for pauses or itnerruptions of the 
    % experiment if they have markers or you know their indices in the data
    for i = 1:numel(start_beep)-1
        removeindices(i+1, :) = [EEG.event(6*i-1).latency+50 EEG.event(6*i+1).latency-50];
    end
    removeindices(end, :) = [EEG.event(end-1).latency+50 EEG.pnts]; % remove from last event to the end
    

    %%
    % filter for plot
    EEG_plot = pop_eegfiltnew(EEG, 'locutoff',0.5, 'hicutoff', 40,'plotfreqz',0);
    
    % plot
    fig1 = figure; set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    plot(normalize(EEG_plot.data') + [1:10:10*EEG_plot.nbchan], 'color', [78 165 216]/255)
    yticks([])
    
    xlim([0 EEG.pnts])
    ylim([-10 10*EEG_plot.nbchan+10])
    
    hold on
    
    % plot lines for valid times
    
    for i = 1:size(removeindices,1)
        plot([removeindices(i,1) removeindices(i,1)],ylim,'k', 'LineStyle','-')
        plot([removeindices(i,2) removeindices(i,2)],ylim,'k', 'LineStyle','--')
    end
    title(['Subject ', num2str(subject), ', Non-Exp Segments: From Solid Line to Next Dashed Line on the Right'])
    
    % save plot
    print(gcf,fullfile(input_filepath,[bemobil_config.filename_prefix num2str(subject) '_raw-full_EEG.png']),'-dpng')
    close

    %% Add other events like flexion and extension
    % Important Note: before rejecting the non-exp segments you must add
    % the flexion/Extension start events

    % cd 'C:\Morteza\Analysis\ANSYMB2024\data\Trials_info\sub-7'
    % load('subj_7_Trials_encoder_events.mat')
    % cd 'C:\Morteza\Analysis\ANSYMB2024\Code\data_processing'

    

    %% reject
    EEG = eeg_eegrej(EEG, removeindices);   
    
    %% processing wrappers for basic processing and AMICA
    
    % do basic preprocessing, line noise removal, and channel interpolation
	[ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_process_all_EEG_preprocessing(subject, bemobil_config, ALLEEG, EEG, CURRENTSET, force_recompute);

    % start the processing pipeline for AMICA
	bemobil_process_all_AMICA(ALLEEG, EEG_preprocessed, CURRENTSET, subject, bemobil_config, force_recompute);

end

bemobil_copy_plots_in_one(bemobil_config)

subjects
subject

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')


