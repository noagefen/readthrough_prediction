clear;
load('Y_rand.mat');

% Configuration
genes_per_file = 10;

% Sizes
num_of_genes = size(Y_rand, 1);
rand_num     = size(Y_rand, 3);

% Output directory
outdir = 'rand_stopwin_fastas';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

for n = 1:rand_num
    % Number of output files for this slice
    num_files = ceil(num_of_genes / genes_per_file);

    for f = 1:num_files
        start_idx = (f-1)*genes_per_file + 1;
        end_idx   = min(f*genes_per_file, num_of_genes);
        iter_num  = end_idx - start_idx + 1;

        % Preallocate struct array with required fields for fastawrite
        Y_stop_win = struct('Sequence', cell(iter_num,1), 'Header', cell(iter_num,1));

        for j = 1:iter_num
            idx = start_idx + j - 1;

            % Robust header extraction (handles cell-of-cell vs. direct char/string)
            hdrCell = Y_rand{idx, 1, n};
            if iscell(hdrCell)
                header = hdrCell{1};
            else
                header = hdrCell;
            end

            % ORF and UTR3 for this slice
            orf  = upper(char(Y_rand{idx, 2, n}));
            utr3 = upper(char(Y_rand{idx, 3, n}));

            % Take last 150 nt of ORF (or full if shorter)
            if length(orf) >= 150
                orf_150 = orf(end-149:end);
            else
                orf_150 = orf;
            end

            % Take first 150 nt of UTR3 (or full if shorter)
            if length(utr3) >= 150
                utr3_150 = utr3(1:150);
            else
                utr3_150 = utr3;
            end

            % Build FASTA entry
            Y_stop_win(j).Sequence = [orf_150 utr3_150]; % concatenated char array
            Y_stop_win(j).Header   = char(header);       % ensure char for fastawrite
        end

        % File name: Yrand<n>_stopwin_<start>_<end>.fas
        fname = sprintf('Yrand%d_stopwin_%d_%d.fas', n, start_idx, end_idx);
        fullpath = fullfile(outdir, fname);

        % Write FASTA
        fastawrite(fullpath, Y_stop_win);
    end
end