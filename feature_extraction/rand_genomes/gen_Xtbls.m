clear;

%% -------------------- FILTER INDEX --------------------
load('../wt_25C_tbl_orig.mat');

fileID = fopen('../reference_set_2693_genes.txt');
filtered_genes = textscan(fileID, '%s');
fclose(fileID);

% wt_25C
wt_25C_index = [];
for i = 1:height(wt_25C_tbl_orig)
     if ~isempty(find([filtered_genes{:}] == wt_25C_tbl_orig.gene(i),1)) && ...
             wt_25C_tbl_orig.cds_RPKM(i)>5 && wt_25C_tbl_orig.ext_RPKM(i)>0.5
         wt_25C_index = [wt_25C_index;i]; %#ok<AGROW>
    end
end
keep = wt_25C_index(:);  % column vector

%% -------------------- CONFIG --------------------
num_slices = 20;
base_csv   = '../WT25Data/Xtbl_wt25.csv';

% Helper for loading vars with expected names:
getVar = @(S, nm) pickVarOrOnly(S, nm);

% Read base once
base_all = readtable(base_csv);

% If base is NOT already filtered to the 2693 rows, and it’s the full set,
% you can uncomment the next line to filter the base too.
% base_all = base_all(keep, :);

%% -------------------- MAIN LOOP --------------------
for num = 1:num_slices
    T = base_all;  % fresh copy of the base for this dataset

    %% --- CAI: CAI_table_<num>_reorg.mat -> replace CAI_30/50/70 ---
    S = tryload(sprintf('CAI_table_%d_reorg.mat', num));
    if ~isempty(S)
        Tsrc = getVar(S, sprintf('CAI_table_%d_reorg', num));
        Tsrc = filterToKeep(Tsrc, keep, height(T), sprintf('CAI_table_%d_reorg', num));
        T = assignCols(T, Tsrc, {'CAI_30','CAI_50','CAI_70'}, {'CAI_30','CAI_50','CAI_70'});
    end

    %% --- tAI: tAI_table_<num>_reorg.mat -> replace tAI_30/50/70 ---
    S = tryload(sprintf('tAI_table_%d_reorg.mat', num));
    if ~isempty(S)
        Tsrc = getVar(S, sprintf('tAI_table_%d_reorg', num));
        Tsrc = filterToKeep(Tsrc, keep, height(T), sprintf('tAI_table_%d_reorg', num));
        T = assignCols(T, Tsrc, {'tAI_30','tAI_50','tAI_70'}, {'tAI_30','tAI_50','tAI_70'});
    end

    %% --- orig3utrs MFE: orig3utrs_X_mfe_tbl_<num>_reorg.mat -> MFE, dist_bp ---
    S = tryload(sprintf('orig3utrs_X_mfe_tbl_%d_reorg.mat', num));
    if ~isempty(S)
        Tsrc = getVar(S, sprintf('orig3utrs_X_mfe_tbl_%d_reorg', num));
        Tsrc = filterToKeep(Tsrc, keep, height(T), sprintf('orig3utrs_X_mfe_tbl_%d_reorg', num));
        T = assignCols(T, Tsrc, {'mfe','mfe_indx'}, {'MFE','dist_bp'});
    end

    %% --- stop-window MFE: stop_win_X_mfe_tbl_<num>_reorg.mat -> stop_win_mfe, stop_win_mfe_indx ---
    S = tryload(sprintf('stop_win_X_mfe_tbl_%d_reorg.mat', num));
    if ~isempty(S)
        Tsrc = getVar(S, sprintf('stop_win_X_mfe_tbl_%d_reorg', num));
        Tsrc = filterToKeep(Tsrc, keep, height(T), sprintf('stop_win_X_mfe_tbl_%d_reorg', num));
        T = assignCols(T, Tsrc, {'stop_win_mfe','stop_win_mfe_indx'}, {'stop_win_mfe','stop_win_mfe_indx'});
    end

    %% --- X_utr3stop: X_utr3stop_tbl_<num>.mat -> nis_stop ---
    S = tryload(sprintf('X_utr3stop_tbl_%d_reorg.mat', num));
    if ~isempty(S)
        Tsrc = getVar(S, sprintf('X_utr3stop_tbl_%d_reorg', num));
        Tsrc = filterToKeep(Tsrc, keep, height(T), sprintf('X_utr3stop_tbl_%d_reorg', num));
        T = assignCols(T, Tsrc, {'X_utr3stop'}, {'nis_stop'});
    end

    %% --- Codon features: codon_features_<num>.mat -> codon_m_6...codon_m_1 ---
    S = tryload(sprintf('codon_features_%d_reorg.mat', num));
    if ~isempty(S)
        Tsrc = getVar(S, sprintf('codon_features_%d_reorg', num));
        Tsrc = filterToKeep(Tsrc, keep, height(T), sprintf('codon_features_%d_reorg', num));
        % Assign if present in base
        codonCols = {'codon_m_6','codon_m_5','codon_m_4','codon_m_3','codon_m_2','codon_m_1'};
        T = assignCols(T, Tsrc, codonCols, codonCols);
    end

    %% --- Y_features (nt_m01..nt_m18 from nt_m_1..nt_m_18; nt_p04..nt_p16 from nt_p_4..nt_p_16) ---
    S = tryload(sprintf('Y_features_%d_reorg.mat', num));
    if ~isempty(S)
        Tsrc = getVar(S, sprintf('Y_features_%d_reorg', num));
        Tsrc = filterToKeep(Tsrc, keep, height(T), sprintf('Y_features_%d_reorg', num));

        % nt_m: map nt_m_1..18 -> nt_m01..18 (zero-padded for 1..9)
        for k = 1:18
            src = sprintf('nt_m_%d', k);
            if k < 10
                dst = sprintf('nt_m0%d', k);
            else
                dst = sprintf('nt_m%d', k);
            end
            if hasVar(T, dst) && hasVar(Tsrc, src)
                T.(dst) = Tsrc.(src);
            else
                warnMissingColumn(dst, src, 'Y_features', num);
            end
        end

        % nt_p: map nt_p_4..16 -> nt_p04..16
        for k = 4:16
            src = sprintf('nt_p_%d', k);
            if k < 10
                dst = sprintf('nt_p0%d', k);
            else
                dst = sprintf('nt_p%d', k);
            end
            if hasVar(T, dst) && hasVar(Tsrc, src)
                T.(dst) = Tsrc.(src);
            else
                warnMissingColumn(dst, src, 'Y_features', num);
            end
        end
    end

    %% -------- Save per-dataset combined table --------
    outVar = sprintf('Xtbl_rand_%d', num);
    Out = struct();
    Out.(outVar) = T;
    save(sprintf('%s.mat', outVar), '-struct', 'Out');
end

%% ==================== HELPERS ====================
function S = tryload(path)
    if isfile(path)
        S = load(path);
    else
        warning('Missing file: %s', path);
        S = [];
    end
end

function T = pickVarOrOnly(S, expectedName)
    if isfield(S, expectedName)
        T = S.(expectedName);
        return;
    end
    fns = fieldnames(S);
    if numel(fns) == 1
        T = S.(fns{1});
        warning('Variable %s not found; using %s instead.', expectedName, fns{1});
    else
        % Prefer a table if present
        tblIdx = find(structfun(@istable, S), 1);
        if ~isempty(tblIdx)
            fn = fns{tblIdx};
            T = S.(fn);
            warning('Variable %s not found; using table %s instead.', expectedName, fn);
        else
            error('Variable %s not found and multiple variables in MAT.', expectedName);
        end
    end
end

function Tsrc = filterToKeep(Tsrc, keep, targetHeight, label)
    % Filters source table to rows 'keep' iff needed to match target height.
    % If heights already match, returns as-is.
    if height(Tsrc) == targetHeight
        return;
    end
    if max(keep) <= height(Tsrc)
        Tsrc = Tsrc(keep, :);
    else
        warning('%s: cannot filter to keep (max index %d > height %d).', label, max(keep), height(Tsrc));
    end
end

function tf = hasVar(T, name)
    tf = any(strcmp(T.Properties.VariableNames, name));
end

function T = assignCols(T, Tsrc, srcNames, dstNames)
    % Assign columns by name with simple type compatibility.
    for ii = 1:numel(srcNames)
        src = srcNames{ii};
        dst = dstNames{ii};
        if hasVar(Tsrc, src)
            if ~hasVar(T, dst)
                warning('Destination column %s not found in base table.', dst);
                continue;
            end
            vals = Tsrc.(src);
            % Ensure column length matches height(T)
            if numel(vals) ~= height(T)
                warning('Length mismatch for %s -> %s (%d vs %d). Skipping assignment.', src, dst, numel(vals), height(T));
                continue;
            end
            T.(dst) = vals;
        else
            warning('Source column %s not found in source table.', src);
        end
    end
end

function warnMissingColumn(dst, src, block, num)
    warning('%s_%d: %s (dst) or %s (src) not found; skipped.', block, num, dst, src);
end
