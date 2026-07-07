% X table including new feature. only genes with wt25 RE valid data.
% Y table for wt25 valid genes.

clear;
load('wt_25C_tbl_orig.mat');

fileID = fopen('reference_set_2693_genes.txt');
filtered_genes = textscan(fileID, '%s');

% wt_25C
wt_25C_index = [];
for i = 1:height(wt_25C_tbl_orig)
     if ~isempty(find([filtered_genes{:}] == wt_25C_tbl_orig.gene(i),1)) && ...
             wt_25C_tbl_orig.cds_RPKM(i)>5 && wt_25C_tbl_orig.ext_RPKM(i)>0.5
         wt_25C_index = [wt_25C_index;i];
    end
end

features_tbl_orig = readtable('feature_file.csv');
X_tbl_orig = features_tbl_orig(:,[1,152,122:127,131,130,129,128,132:137,141,140,...
    139,138,142:147,151,150,149,148,154:159,163,162,161,160,120,121,71,70,69,...
    68,67,66,65,64,63,62,61,60,59,58,57,56,55,54,6,41:53,5,8,12,13,29,30]);


load('CAI_table_reorg.mat');
load('codon_features_reorg.mat');
load('similar_to_stop_reorg.mat');
load('stop_win_X_mfe_tbl_reorg.mat');
load('tAI_table_reorg.mat');

X_w_new = [X_tbl_orig , CAI_table_reorg(:,2:end) , codon_features_reorg(:,2:end) ,...
    similar_to_stop_reorg(:,2:end) , stop_win_X_mfe_tbl_reorg , tAI_table_reorg(:,2:end)];

Xtbl_wt25 = X_w_new(wt_25C_index,:);
% Ytbl_wt25 = wt_25C_tbl_orig.RE(wt_25C_index);
Ytbl_wt25 = wt_25C_tbl_orig(wt_25C_index,1:2);

save Xtbl_wt25.mat Xtbl_wt25;
save Ytbl_wt25.mat Ytbl_wt25;

writetable(Xtbl_wt25,'WT25Data/Xtbl_wt25.csv');
writetable(Ytbl_wt25,'WT25Data/Ytbl_wt25.csv');