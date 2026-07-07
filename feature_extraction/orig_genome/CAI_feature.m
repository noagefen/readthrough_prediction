clear;
load('CAI_w_map.mat');
load('Y_genes_orig3utrs.mat');

CAI_30 = zeros(height(Y_genes_orig3utrs),1);
CAI_50 = zeros(height(Y_genes_orig3utrs),1);
CAI_70 = zeros(height(Y_genes_orig3utrs),1);
genes = strings(height(Y_genes_orig3utrs),1);

for i = 1:height(Y_genes_orig3utrs)
    genes(i) = Y_genes_orig3utrs{i,1}{1};
    seq = Y_genes_orig3utrs{i,2};
    cai_weights = zeros(1,length(seq)/3-1); %skip start codon
    len_codons = length(cai_weights);
    for j = 1:len_codons
        codon = seq(j*3+1:j*3+3);
        cai_weights(j) = CAI_w_map(codon);
    end

    if len_codons < 30
        CAI_30(i) = geomean(cai_weights);
    else 
        CAI_30(i) = geomean(cai_weights(end-29:end));
    end

    if len_codons < 50
        CAI_50(i) = geomean(cai_weights);
    else 
        CAI_50(i) = geomean(cai_weights(end-49:end));
    end

    if len_codons < 70
        CAI_70(i) = geomean(cai_weights);
    else 
        CAI_70(i) = geomean(cai_weights(end-69:end));
    end
end

CAI_table = [table(genes) table(CAI_30) table(CAI_50) table(CAI_70)];
save('CAI_table.mat', 'CAI_table', '-v7.3');