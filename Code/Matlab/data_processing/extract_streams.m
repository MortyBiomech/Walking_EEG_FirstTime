function out = extract_streams(subject_id, rawdata_path, study_identifier)

    %% Load XDF files (it takes a while to load all the streams!)
    streams    = xdf_load_matlab(subject_id, rawdata_path, study_identifier);
    
    %% Concatenating Signals
    % Identify the indeces for EEG, EMG, and Encoder(Exp_data) signals
    EEG_indx = [];
    EMG_indx = [];
    GRF_indx = [];
    
    for i = 1:length(streams)
        for j = 1:4 % 4: EEG_trigger_markers, EEG, EMG, GRF
            type = streams{1, i}{1, j}.info.type;
            if strcmp(type, 'EEG')
                EEG_indx = cat(2, EEG_indx, j);
            end
            if strcmp(type, 'EMG')
                EMG_indx = cat(2, EMG_indx, j);
            end
            if strcmp(type, 'Force')
                GRF_indx = cat(2, GRF_indx, j);
            end
        end
    end
    
    
    % Initialize cell arrays to store data
    temp_EEG      = cell(1, numel(EEG_indx));
    temp_EEG_time = cell(1, numel(EEG_indx));
    
    temp_EMG      = cell(1, numel(EMG_indx));
    temp_EMG_time = cell(1, numel(EMG_indx));
    
    temp_GRF      = cell(1, numel(GRF_indx));
    temp_GRF_time = cell(1, numel(GRF_indx));
    
    % Loop to collect data into cell arrays
    for i = 1:length(streams)
        temp_EEG{i} = streams{1, i}{1, EEG_indx(i)}.time_series;
        temp_EEG_time{i} = streams{1, i}{1, EEG_indx(i)}.time_stamps;

        temp_EMG{i} = streams{1, i}{1, EMG_indx(i)}.time_series;
        temp_EMG_time{i} = streams{1, i}{1, EMG_indx(i)}.time_stamps;

        temp_GRF{i} = streams{1, i}{1, GRF_indx(i)}.time_series;
        temp_GRF_time{i} = streams{1, i}{1, GRF_indx(i)}.time_stamps;
    end
    

    % Concatenate cell arrays into final arrays
    All_EEG      = double([temp_EEG{:}]);
    All_EEG_time = [temp_EEG_time{:}];
    
    % Check if the EMG sensor ids are the same in all sessions
    for i = 1:length(temp_EMG)
        if size(temp_EMG{1, i},1) ~= size(temp_EMG{1, 1},1)
            % Instruction message
            prompt = {sprintf(['New EMG sensors were utilized at ', ...
                'session ', num2str(i), '. Please write down the new ', ...
                'sensor IDs and their corresponding old ones ' ...
                '(e.x. [12 2; 13 3; 14 4; 15 5], this means sensors ', ...
                '12, 13, 14, and 15 were used instead of sensors ', ...
                '2, 3, 4, and 5, respectively):\n'])};
            dlg_title = 'EMG sensors replacement';
            num_lines = 1;
            defaultans = {'[12 2; 13 3; 14 4; 15 5]'};
            
            % Display input dialog
            try 
                answer = ...
                    inputdlg(prompt, dlg_title, num_lines, defaultans);
                answer = str2double(answer{1,1});
                temp_EMG{1, i}(answer(:,2),:) = ...
                    temp_EMG{1, i}(answer(:,1),:);
                temp_EMG{1, i}(answer(:,1),:) = [];
            catch
                disp('No input from user!')
            end
        end
    end

    All_EMG      = double([temp_EMG{:}]);
    All_EMG_time = [temp_EMG_time{:}];
    
    All_GRF      = double([temp_GRF{:}]);
    All_GRF_time = [temp_GRF_time{:}];


    out = struct('All_EEG', All_EEG, 'All_EEG_time', All_EEG_time, ...
        'All_EMG', All_EMG, 'All_EMG_time', All_EMG_time, ...
        'All_GRF', All_GRF, 'All_GRF_time', All_GRF_time);

end
