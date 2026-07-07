clear;

% ----- CONFIG (match your original) -----
genes_per_file = 10;
num_of_genes   = 6686;
num_slices     = 20;   % Yrand1..Yrand20
data_dir       = 'rand_orig3utrs_mfe_res';
% ----------------------------------------

num_files_per_slice = ceil(num_of_genes / genes_per_file);

for n = 1:num_slices
    % Preallocate [mfe, index] like original (one dataset per slice)
    orig3utrs_X_mfe = zeros(num_of_genes, 2);

    for i = 1:num_files_per_slice
        % Determine how many genes in this chunk (last chunk may be short)
        if i == num_files_per_slice
            iter_num = num_of_genes - (i-1)*genes_per_file;
        else
            iter_num = genes_per_file;
        end

        start_idx = (i-1)*genes_per_file + 1;
        end_idx   = (i-1)*genes_per_file + iter_num;

        % File name pattern: Yrand<n>_utr3_<start>_<end>_res.fas
        file_num = sprintf('%d_%d', start_idx, end_idx);
        filename = fullfile(data_dir, sprintf('Yrand%d_utr3_%s_res.fas', n, file_num));

        if ~isfile(filename)
            warning('Missing file for slice %d: %s', n, filename);
            % Leave corresponding rows as zeros; continue to next chunk
            continue;
        end

        % Read FASTA-like results; second output contains the text payload per record
        [~, mfe_res] = fastaread(filename);

        % Defensive: if file has fewer/more entries than expected, cap to available
        m = min(iter_num, numel(mfe_res));
        if m < iter_num
            warning('File %s has %d records; expected %d. Using %d.', filename, numel(mfe_res), iter_num, m);
        end

        % Parse each record: "max energy = %s index = %s"
        for j = 1:m
            % mfe_res{j} is a single line/string like: "max energy = <val> index = <val>"
            tokens = textscan(mfe_res{j}, 'max energy = %s index = %s');
            % Convert to double and place into the correct absolute row
            row_idx = (i-1)*genes_per_file + j;
            orig3utrs_X_mfe(row_idx, 1) = str2double(tokens{1}{1});  % mfe
            orig3utrs_X_mfe(row_idx, 2) = str2double(tokens{2}{1});  % index
        end
    end

    % Build table (same column names as original)
    orig3utrs_X_mfe_tbl = array2table(orig3utrs_X_mfe, ...
        'VariableNames', {'mfe','mfe_indx'});

    % Save with variable name matching the file (no eval)
    varname = sprintf('orig3utrs_X_mfe_tbl_%d', n);
    S = struct();
    S.(varname) = orig3utrs_X_mfe_tbl;
    save([varname '.mat'], '-struct', 'S');
end
