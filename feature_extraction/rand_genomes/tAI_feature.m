clear;
load('tAI_w_map.mat');
load('Y_rand.mat');

dims = size(Y_rand);
rand_num = dims(3);
num_of_genes = size(Y_rand,1);

for n = 1:rand_num
    tAI_30 = zeros(num_of_genes,1);
    tAI_50 = zeros(num_of_genes,1);
    tAI_70 = zeros(num_of_genes,1);
    genes = strings(num_of_genes,1);

    for i = 1:num_of_genes
        genes(i) = Y_rand{i,1,n}{1};
        seq = Y_rand{i,2,n};
        codon_count = floor(length(seq)/3) - 1; % skip start codon
        tai_weights = zeros(1, codon_count);

        for j = 1:codon_count
            codon = seq(j*3+1:j*3+3);
            tai_weights(j) = tAI_w_map(codon);
        end

        % Compute tAI for windows
        tAI_30(i) = geomean(tai_weights(max(1, end-29):end));
        tAI_50(i) = geomean(tai_weights(max(1, end-49):end));
        tAI_70(i) = geomean(tai_weights(max(1, end-69):end));
    end

    % Build and save table
    tAI_table = table(genes, tAI_30, tAI_50, tAI_70);
    varname = sprintf('tAI_table_%d', n);
    eval([char(varname) ' = tAI_table;']);
    save([char(varname) '.mat'], char(varname), '-v7.3');
end
