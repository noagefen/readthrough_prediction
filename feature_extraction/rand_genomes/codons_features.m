clear;
load('Y_rand.mat');

codons_m_num = 6;
dims = size(Y_rand);
rand_num = dims(3);
num_of_genes = size(Y_rand,1);

for n = 1:rand_num
    gene = Y_rand(:,1,n);
    codons = strings(num_of_genes, codons_m_num);

    for i = 1:num_of_genes
        orf = upper(Y_rand{i,2,n});
        for j = 1:codons_m_num
            minus_start = j*3+2;
            minus_end = minus_start-2;
            if length(orf) >= minus_start
                codons(i,j) = orf(end-minus_start:end-minus_end);
            end
        end
    end

    % Build table
    codon_features = table(gene);
    for j = codons_m_num:-1:1
        colname = sprintf("codon_m_%d", j);
        codon_features.(colname) = codons(:,j);
    end

    % Assign dynamic variable name and save
    varname = sprintf('codon_features_%d', n);
    eval([varname ' = codon_features;']);
    save([varname '.mat'], varname);
end
