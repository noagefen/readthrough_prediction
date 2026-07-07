clear;
load('Y_genes.mat');

codons_m_num = 6;
for j = 1:codons_m_num
    eval(['codon_m_' num2str(j) ' = strings(height(Y_genes),1);']);
end

gene = Y_genes(:,1);

for i = 1:height(Y_genes)
    orf = upper(Y_genes{i,2});
    for j = 1:codons_m_num
        minus_start = j*3+2; minus_end = minus_start-2;
        eval(['codon_m_' num2str(j) '(i) = orf(end-minus_start:end-minus_end);']);
    end
end

codon_features = table(gene);
for j = codons_m_num:-1:1
    eval(['codon_features = [codon_features array2table(codon_m_' num2str(j) ')];']);
end
save codon_features.mat codon_features;