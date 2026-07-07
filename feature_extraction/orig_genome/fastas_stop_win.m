clear;
load('Y_genes_orig3utrs.mat');

genes_per_file = 10;
num_of_genes = 6686;

for i = 1:ceil(num_of_genes/genes_per_file)
    if i==ceil(num_of_genes/genes_per_file)
        iter_num = num_of_genes-(i-1)*genes_per_file;
    else
        iter_num = genes_per_file;
    end
    Y_stop_win = struct('Sequence', cell(iter_num,1), 'Header', cell(iter_num,1));
    for j = 1:iter_num
        utr3 = Y_genes_orig3utrs{(i-1)*genes_per_file+j,3};
        orf = Y_genes_orig3utrs{(i-1)*genes_per_file+j,2};
        if length(orf)>=150
            orf_150 = orf(end-149:end);
        else
            orf_150 = orf;
        end
        if length(utr3)>=150
            utr3_150 = utr3(1:150);
        else
            utr3_150 = utr3;
        end
        Y_stop_win(j).Sequence = strcat(orf_150,utr3_150);
        Y_stop_win(j).Header = Y_genes_orig3utrs{(i-1)*genes_per_file+j,1}{1};
    end
    filename = strcat('stop_win_fastas/Y_stop_win_',num2str((i-1)*genes_per_file+1),'_',...
        num2str((i-1)*genes_per_file+iter_num),'.fas');
    fastawrite(filename,Y_stop_win);
end