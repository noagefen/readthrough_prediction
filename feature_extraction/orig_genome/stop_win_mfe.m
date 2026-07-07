genes_per_file = 10;
num_of_genes = 6686;
stop_win_X_mfe = zeros(num_of_genes,2);

for i = 1:ceil(num_of_genes/genes_per_file)
    if i==ceil(num_of_genes/genes_per_file)
        iter_num = num_of_genes-(i-1)*genes_per_file;
    else
        iter_num = genes_per_file;
    end
    file_num = strcat(num2str((i-1)*genes_per_file+1),'_',...
        num2str((i-1)*genes_per_file+iter_num));
    filename = strcat('stop_win_mfe/Y_stop_win_',file_num,'_res.fas');
    [~,mfe_res] = fastaread(filename);
    for j = 1:iter_num
        mfe = textscan(mfe_res{j},'max energy = %s index = %s');
        stop_win_X_mfe((i-1)*genes_per_file+j,1) = str2double(mfe{1}{1});
        stop_win_X_mfe((i-1)*genes_per_file+j,2) = str2double(mfe{2}{1});
    end
end

stop_win_X_mfe_tbl = array2table(stop_win_X_mfe,'VariableNames',{'stop_win_mfe','stop_win_mfe_indx'});
save stop_win_X_mfe_tbl.mat stop_win_X_mfe_tbl;