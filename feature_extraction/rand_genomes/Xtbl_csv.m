% ============================================
% Convert Xtbl_rand_1.mat ... Xtbl_rand_20.mat 
% into CSV files
% ============================================
clear;

inputFolder = '.';              % folder where your .mat files are
outputFolder = './Xtbl_csv';  % output folder

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

for i = 1:20
    matFile = fullfile(inputFolder, sprintf('Xtbl_rand_%d.mat', i));
    fprintf('Loading %s...\n', matFile);

    % Load .mat file
    S = load(matFile);

    % Identify the actual variable inside the .mat file
    vars = fieldnames(S);
    if length(vars) ~= 1
        warning('Multiple variables found in %s. Using the first one.\n', matFile);
    end
    varName = vars{1};
    data = S.(varName);

    % Convert to table if needed
    if istable(data)
        T = data;
    elseif isstruct(data)
        % struct array or a MATLAB table saved as struct
        try
            T = struct2table(data);
        catch
            error('Cannot convert struct in %s to a table.', matFile);
        end
    elseif isnumeric(data)
        % numeric matrix → table with generic column names
        colNames = arrayfun(@(k) sprintf('col_%d', k), 1:size(data, 2), 'UniformOutput', false);
        T = array2table(data, 'VariableNames', colNames);
    else
        error('Unsupported data type in %s.', matFile);
    end

    % Save to CSV
    outFile = fullfile(outputFolder, sprintf('Xtbl_rand_%d.csv', i));
    writetable(T, outFile);

    fprintf('Saved → %s\n', outFile);
end

fprintf('\n🎉 All .mat files converted successfully!\n');
    