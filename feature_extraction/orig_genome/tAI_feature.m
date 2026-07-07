clear;
load('tAI_w_map.mat');
load('Y_genes_orig3utrs.mat');

tAI_30 = zeros(height(Y_genes_orig3utrs),1);
tAI_50 = zeros(height(Y_genes_orig3utrs),1);
tAI_70 = zeros(height(Y_genes_orig3utrs),1);
genes = strings(height(Y_genes_orig3utrs),1);

for i = 1:height(Y_genes_orig3utrs)
    genes(i) = Y_genes_orig3utrs{i,1}{1};
    seq = Y_genes_orig3utrs{i,2};
    tai_weights = zeros(1,length(seq)/3-1); %skip start codon
    len_codons = length(tai_weights);
    for j = 1:len_codons
        codon = seq(j*3+1:j*3+3);
        tai_weights(j) = tAI_w_map(codon);
    end

    if len_codons < 30
        tAI_30(i) = geomean(tai_weights);
    else 
        tAI_30(i) = geomean(tai_weights(end-29:end));
    end

    if len_codons < 50
        tAI_50(i) = geomean(tai_weights);
    else 
        tAI_50(i) = geomean(tai_weights(end-49:end));
    end

    if len_codons < 70
        tAI_70(i) = geomean(tai_weights);
    else 
        tAI_70(i) = geomean(tai_weights(end-69:end));
    end
end

tAI_table = [table(genes) table(tAI_30) table(tAI_50) table(tAI_70)];
save('tAI_table.mat', 'tAI_table', '-v7.3');