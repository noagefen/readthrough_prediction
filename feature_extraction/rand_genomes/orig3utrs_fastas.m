clear;
load('Y_rand.mat');

% Configuration
genes_per_file = 10;

% Sizes
num_of_genes = 6686;
rand_num     = size(Y_rand, 3);

% Output directory (already chosen)
outdir = 'rand_orig3utrs_fasta';
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

        % Preallocate struct array for fastawrite
        Y_utr3 = struct('Sequence', cell(iter_num,1), 'Header', cell(iter_num,1));

        for j = 1:iter_num
            idx = start_idx + j - 1;

            % Header (robust to cell-of-cell or direct char/string)
            hdrVal = Y_rand{idx, 1, n};
            if iscell(hdrVal)
                header = hdrVal{1};
            else
                header = hdrVal;
            end

            % UTR3 value (robust to string/char/cell)
            utrVal = Y_rand{idx, 3, n};
            if iscell(utrVal)
                utrVal = utrVal{1};
            end

            % Ensure char type for fastawrite
            Y_utr3(j).Sequence = char(utrVal);
            Y_utr3(j).Header   = char(header);
        end

        % File name: Yrand<n>_utr3_<start>_<end>.fas
        fname    = sprintf('Yrand%d_utr3_%d_%d.fas', n, start_idx, end_idx);
        fullpath = fullfile(outdir, fname);

        % Write FASTA (Bioinformatics Toolbox)
        fastawrite(fullpath, Y_utr3);
    end
end
