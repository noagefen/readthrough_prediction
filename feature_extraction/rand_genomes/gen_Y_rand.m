clear;
load ('Y_genes_orig3utrs.mat');
rand_num = 20;
Y_rand = cell(height(Y_genes_orig3utrs),3,rand_num);
for i = 1:height(Y_genes_orig3utrs)
    name = Y_genes_orig3utrs{i,1};
    seq = Y_genes_orig3utrs{i,2};
    utr3 = Y_genes_orig3utrs{i,3};
    for j = 1:rand_num
        Y_rand(i,1,j) = {name};
        Y_rand(i,3,j) = {perm_nt(utr3)};
    end
    rand_seq = perm_codons(seq,rand_num);
    for j = 1:rand_num
        Y_rand(i,2,j) = rand_seq(j);
    end
end
save Y_rand.mat Y_rand -v7.3;