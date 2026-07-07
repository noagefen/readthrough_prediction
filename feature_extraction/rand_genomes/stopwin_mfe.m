clear;

% ---------- CONFIG ----------
data_dir       = 'rand_stopwin_mfe';  % where Yrand*_stopwin_*_res.fas live
genes_per_file = 10;
num_genes      = 6686;                % as in your dataset
num_slices     = 20;                  % Yrand1..Yrand20
% --------------------------------

num_files_per_slice = ceil(num_genes / genes_per_file);

for n = 1:num_slices
    % Preallocate like your original: [mfe, index]
    stop_win_X_mfe = zeros(num_genes, 2);

    for i = 1:num_files_per_slice
        % Determine chunk bounds and size (last chunk is shorter)
        if i == num_files_per_slice
            iter_num = num_genes - (i-1)*genes_per_file;  % 6
        else
            iter_num = genes_per_file;                     % 10
        end
        start_idx = (i-1)*genes_per_file + 1;             % e.g., 6681 for last chunk
        end_idx   = (i-1)*genes_per_file + iter_num;      % e.g., 6686 for last chunk

        % Expected filename: Yrand<n>_stopwin_<start>_<end>_res.fas
        file_num = sprintf('%d_%d', start_idx, end_idx);
        fpath    = fullfile(data_dir, sprintf('Yrand%d_stopwin_%s_res.fas', n, file_num));

        if isfile(fpath)
            % Normal path (unchanged)
            [~, mfe_res] = fastaread(fpath);
            if ischar(mfe_res), mfe_res = {mfe_res}; end
            m = min(iter_num, numel(mfe_res));
            for j = 1:m
                tokens = textscan(mfe_res{j}, 'max energy = %s index = %s');
                stop_win_X_mfe((i-1)*genes_per_file + j, 1) = str2double(tokens{1}{1});
                stop_win_X_mfe((i-1)*genes_per_file + j, 2) = str2double(tokens{2}{1});
            end

        elseif i == num_files_per_slice
            % ----- Specific fix: use oversized final file 6681_6690 -----
            alt_name = sprintf('Yrand%d_stopwin_6681_6690_res.fas', n);
            alt_path = fullfile(data_dir, alt_name);
            if isfile(alt_path)
                [~, mfe_res] = fastaread(alt_path);
                if ischar(mfe_res), mfe_res = {mfe_res}; end

                % Take first 6 records -> rows 6681..6686
                M = min(iter_num, numel(mfe_res));  % iter_num is 6
                for j = 1:M
                    tokens = textscan(mfe_res{j}, 'max energy = %s index = %s');
                    row_idx = start_idx + j - 1;    % 6681..6686
                    if row_idx >= 1 && row_idx <= num_genes
                        stop_win_X_mfe(row_idx, 1) = str2double(tokens{1}{1});
                        stop_win_X_mfe(row_idx, 2) = str2double(tokens{2}{1});
                    end
                end
            end
            % If the alt file is also missing, rows stay zeros (original behavior).
        end
    end

    % Build table (same columns as original)
    stop_win_X_mfe_tbl = array2table(stop_win_X_mfe, ...
        'VariableNames', {'stop_win_mfe','stop_win_mfe_indx'});

    % Save with variable name matching the slice (no eval)
    varname = sprintf('stop_win_X_mfe_tbl_%d', n);
    S = struct();
    S.(varname) = stop_win_X_mfe_tbl;
    save([varname '.mat'], '-struct', 'S');
end
