clear;

% --- Build the new index order (your code) ---
load('../Y_genes_orig3utrs.mat');
features_tbl_orig = readtable('../feature_file.csv');
idxs = zeros(height(features_tbl_orig),1);

for i = 1:height(features_tbl_orig)
    gene = erase(features_tbl_orig.transcript{i},'_mRNA');
    idxs(i) = find(string(Y_genes_orig3utrs(:,1)) == gene, 1, 'first');
end

% (Optional sanity check)
if any(idxs == 0)
    error('Some transcripts from feature_file.csv were not found in Y_genes_orig3utrs.');
end

% ===========================================================
num_slices    = 20;            % dataset numbers: 1..20
in_dir        = '.';           % folder where <feature>_<n>.mat live
out_dir       = '.';           % where to write <feature>_<n>_reorg.mat
feature_list  = { ...
    'CAI_table', ...
    'codon_features', ...
    'orig3utrs_X_mfe_tbl', ...
    'stop_win_X_mfe_tbl', ...
    'tAI_table', ...
    'X_utr3stop_tbl', ...
    'Y_features', ...
};
% ===========================================================

for f = 1:numel(feature_list)
    base = feature_list{f};  % e.g., 'stop_win_X_mfe_tbl'
    for n = 1:num_slices
        in_mat  = fullfile(in_dir,  sprintf('%s_%d.mat', base, n));
        out_mat = fullfile(out_dir, sprintf('%s_%d_reorg.mat', base, n));
        var_in  = sprintf('%s_%d', base, n);
        var_out = sprintf('%s_%d_reorg', base, n);

        if ~isfile(in_mat)
            warning('Missing file: %s (skipping)', in_mat);
            continue;
        end

        % Load input .mat
        S = load(in_mat);

        % Pick the table variable
        if isfield(S, var_in)
            T = S.(var_in);
        else
            % Fallback: if the MAT has a single variable, use it
            fns = fieldnames(S);
            if numel(fns) == 1
                T = S.(fns{1});
                warning('Variable %s not found in %s; used %s instead.', var_in, in_mat, fns{1});
            else
                warning('Variable %s not found in %s and multiple variables exist. Skipping.', var_in, in_mat);
                continue;
            end
        end

        % Ensure we have a table
        if ~istable(T)
            % Try to convert common types to table
            if isstruct(T)
                T = struct2table(T);
            else
                T = array2table(T);  % columns will be Var1,Var2,...
            end
        end

        % Bounds check
        if max(idxs) > height(T)
            warning('Index out of bounds for %s: max(idxs)=%d > height(T)=%d. Skipping.', in_mat, max(idxs), height(T));
            continue;
        end

        % Reorder
        T_reorg = T(idxs, :);

        % Save with exact variable name (no eval)
        Out = struct();
        Out.(var_out) = T_reorg;
        save(out_mat, '-struct', 'Out');
    end
end
