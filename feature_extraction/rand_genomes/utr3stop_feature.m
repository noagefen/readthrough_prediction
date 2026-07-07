clear;
load('Y_rand.mat');

dims = size(Y_rand);
rand_num = dims(3);
num_of_genes = size(Y_rand,1);

for n = 1:rand_num
    X_utr3stop = strings(num_of_genes,1);
    indx = zeros(num_of_genes,1);

    for i = 1:num_of_genes
        utr3 = upper(Y_rand{i,3,n});
        X_utr3stop(i) = "NON";

        for j = 1:3:(length(utr3)-2)
            codon = utr3(j:j+2);
            if codon == "TAA" || codon == "TAG" || codon == "TGA"
                X_utr3stop(i) = codon;
                indx(i) = j;
                break;
            end
        end
    end

    % Convert to table
    X_utr3stop_tbl = array2table(X_utr3stop);

    % Save with variable name matching slice number
    varname = sprintf('X_utr3stop_tbl_%d', n);
    eval([char(varname) ' = X_utr3stop_tbl;']);
    save([char(varname) '.mat'], char(varname));
end
