clear;
load('Y_rand.mat');

dims = size(Y_rand);
rand_num = dims(3);
num_of_genes = size(Y_rand,1);

for n = 1:rand_num
    similar_to_stop_cnt = zeros(num_of_genes,1);
    similar_cnt_30 = zeros(num_of_genes,1);
    similar_cnt_50 = zeros(num_of_genes,1);
    similar_cnt_70 = zeros(num_of_genes,1);
    gene = strings(num_of_genes,1);

    for i = 1:num_of_genes
        gene(i) = Y_rand{i,1,n};
        orf = upper(Y_rand{i,2,n});
        orf_len = length(orf);
        if orf_len < 3, continue; end

        stop_codon = orf(end-2:end);
        prefix = stop_codon(1:2);

        % Full ORF
        for j = 1:3:(orf_len-2)
            if orf(j:j+1) == prefix
                similar_to_stop_cnt(i) = similar_to_stop_cnt(i) + 1;
            end
        end

        % Windows
        if orf_len >= 90
            start30 = orf_len - 89;
        else
            start30 = 1;
        end
        for j = start30:3:(orf_len-2)
            if orf(j:j+1) == prefix
                similar_cnt_30(i) = similar_cnt_30(i) + 1;
            end
        end

        if orf_len >= 150
            start50 = orf_len - 149;
        else
            start50 = 1;
        end
        for j = start50:3:(orf_len-2)
            if orf(j:j+1) == prefix
                similar_cnt_50(i) = similar_cnt_50(i) + 1;
            end
        end

        if orf_len >= 210
            start70 = orf_len - 209;
        else
            start70 = 1;
        end
        for j = start70:3:(orf_len-2)
            if orf(j:j+1) == prefix
                similar_cnt_70(i) = similar_cnt_70(i) + 1;
            end
        end
    end

    % Save result table
    similar_to_stop_tbl = table(gene, similar_to_stop_cnt, ...
        similar_cnt_70, similar_cnt_50, similar_cnt_30);
    varname = sprintf('similar_to_stop_tbl_%d', n);
    eval([char(varname) ' = similar_to_stop_tbl;']);
    save([char(varname) '.mat'], char(varname));
end
