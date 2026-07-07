clear;
load('CAI_w_map.mat');
load('Y_rand.mat');

dims = size(Y_rand);
rand_num = dims(3);
num_of_genes = size(Y_rand,1);

for n = 1:rand_num
    CAI_30 = zeros(num_of_genes,1);
    CAI_50 = zeros(num_of_genes,1);
    CAI_70 = zeros(num_of_genes,1);
    genes = strings(num_of_genes,1);

    for i = 1:num_of_genes
        genes(i) = Y_rand{i,1,n}{1};
        seq = Y_rand{i,2,n};
        codon_count = floor(length(seq)/3) - 1;  % skip start codon
        cai_weights = zeros(1, codon_count);

        for j = 1:codon_count
            codon = seq(j*3+1:j*3+3);
            cai_weights(j) = CAI_w_map(codon);
        end

        % Compute CAI windows
        CAI_30(i) = geomean(cai_weights(max(1, end-29):end));
        CAI_50(i) = geomean(cai_weights(max(1, end-49):end));
        CAI_70(i) = geomean(cai_weights(max(1, end-69):end));
    end

    % Build and save table
    CAI_table = table(genes, CAI_30, CAI_50, CAI_70);
    varname = sprintf("CAI_table_%d", n);
    eval([char(varname) ' = CAI_table;']);
    save([char(varname) '.mat'], char(varname), '-v7.3');
end