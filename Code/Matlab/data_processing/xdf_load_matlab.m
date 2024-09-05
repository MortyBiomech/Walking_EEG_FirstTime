function streams = xdf_load_matlab(subject_id, rawdata_path, study_identifier)
    
    files_path = [rawdata_path, 'sub-', num2str(subject_id), filesep, ...
        'ses-S00', num2str(study_identifier), filesep, 'eeg', filesep];

    %% Find the numbers of sessions and load streams
    % Get the list of all items in the folder
    items = dir(files_path);
    streams = repmat(cell(1,1), 1, length(items));
   
    % Loop through the items and load xdf files
    for k = 1:length(items)
        
        % Check if the item is a directory and not '.' or '..' and does not
        % contain 'old' in its name.
        if ~strcmp(items(k).name, '.') && ~strcmp(items(k).name, '..') ...
                && ~contains(items(k).name, 'old')

            streams{1,k} = load_xdf(fullfile(files_path, items(k).name));
 
        end

    end

    %% Remove empty cell members
    % Create a logical array where true indicates an empty cell
    emptyCells = cellfun(@isempty, streams);
    
    % Remove the empty cells
    streams(emptyCells) = [];

 
end