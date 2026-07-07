load('Y_genes_orig3utrs.mat');

similar_to_stop_cnt = zeros(height(Y_genes_orig3utrs),1);
similar_cnt_30 = zeros(height(Y_genes_orig3utrs),1);
similar_cnt_50 = zeros(height(Y_genes_orig3utrs),1);
similar_cnt_70 = zeros(height(Y_genes_orig3utrs),1);
gene = strings(height(Y_genes_orig3utrs),1);

for i = 1:height(Y_genes_orig3utrs)
    gene(i) = Y_genes_orig3utrs{i,1};
    orf = upper(Y_genes_orig3utrs{i,2});
    stop_codon = orf(end-2:end);
    for j = 1:3:length(orf)
        if orf(j:j+1) == stop_codon(1:2)
            similar_to_stop_cnt(i) = similar_to_stop_cnt(i) + 1;
        end
    end
    if length(orf) < 90
        similar_cnt_30(i) = similar_to_stop_cnt(i);
    else
        for j = length(orf)-89:3:length(orf)
            if orf(j:j+1) == stop_codon(1:2)
                similar_cnt_30(i) = similar_cnt_30(i) + 1;
            end
        end
    end
    if length(orf) < 150
        similar_cnt_50(i) = similar_to_stop_cnt(i);
    else
        for j = length(orf)-149:3:length(orf)
            if orf(j:j+1) == stop_codon(1:2)
                similar_cnt_50(i) = similar_cnt_50(i) + 1;
            end
        end
    end
    if length(orf) < 210
        similar_cnt_70(i) = similar_to_stop_cnt(i);
    else
        for j = length(orf)-209:3:length(orf)
            if orf(j:j+1) == stop_codon(1:2)
                similar_cnt_70(i) = similar_cnt_70(i) + 1;
            end
        end
    end
end

similar_to_stop_tbl = table(gene,similar_to_stop_cnt,similar_cnt_70,similar_cnt_50,similar_cnt_30);
save similar_to_stop.mat similar_to_stop_tbl;