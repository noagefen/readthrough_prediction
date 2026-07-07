clear
load('Y_rand.mat')
dims = size(Y_rand);
rand_num = dims(3);
nt_m_num = 18; 
nt_p_num = 13;

% Initialize a cell array to store all Y_features tables
%Y_features_all = cell(rand_num, 1);

for n = 1:rand_num
    gene = Y_rand(:,1,n);
    nt_m = strings(height(Y_rand), nt_m_num);
    nt_p = strings(height(Y_rand), nt_p_num);
    stop_codon = strings(height(Y_rand),1);
    utr3_len = zeros(height(Y_rand),1);
    p_site_aa = strings(height(Y_rand),1);
    e_site_aa = strings(height(Y_rand),1);

    for i = 1:height(Y_rand)
        orf = upper(Y_rand{i,2,n});
        utr3 = upper(Y_rand{i,3,n});
        utr3_len(i) = length(utr3);
        e_site_aa(i) = nt2aa(orf(end-8:end-6));
        p_site_aa(i) = nt2aa(orf(end-5:end-3));
        for j = 1:nt_m_num
            nt_m(i,j) = orf(end-j-2);
        end
        for j = 1:min(nt_p_num, utr3_len(i))
            nt_p(i,j) = utr3(j);
        end
        stop_codon(i) = orf(end-2:end);
    end

    % Create table
    Y_features = table(gene,e_site_aa,p_site_aa);
    for j = nt_m_num:-1:1
        colname = sprintf("nt_m_%d", j);
        Y_features.(colname) = nt_m(:,j);
    end
    Y_features.stop_codon = stop_codon;
    for j = 1:nt_p_num
        colname = sprintf("nt_p_%d", j+3);
        Y_features.(colname) = nt_p(:,j);
    end
    Y_features.utr3_len = utr3_len;

    % Create a dynamic variable name and assign it
    varname = sprintf('Y_features_%d', n);
    eval([varname ' = Y_features;']);  % assigns the table to Y_features_n

    % Save the variable with the same name
    save(sprintf('%s.mat', varname), varname);

end
