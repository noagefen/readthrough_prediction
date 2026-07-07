clear;

load('Y_genes_orig3utrs.mat');

load('CAI_table.mat');
load('codon_features.mat');
load('similar_to_stop.mat');
load('stop_win_X_mfe_tbl.mat');
load('tAI_table.mat');

features_tbl_orig = readtable('feature_file.csv');
idxs = zeros(height(features_tbl_orig),1);

for i = 1:height(features_tbl_orig)
    gene = erase(features_tbl_orig.transcript{i},'_mRNA');
    idxs(i) = find(string(Y_genes_orig3utrs(:,1)) == gene);
end

CAI_table_reorg = CAI_table(idxs,:);
codon_features_reorg = codon_features(idxs,:);
similar_to_stop_reorg = similar_to_stop_tbl(idxs,:);
stop_win_X_mfe_tbl_reorg = stop_win_X_mfe_tbl(idxs,:);
tAI_table_reorg = tAI_table(idxs,:);

save CAI_table_reorg.mat CAI_table_reorg;
save codon_features_reorg.mat codon_features_reorg;
save similar_to_stop_reorg.mat similar_to_stop_reorg;
save stop_win_X_mfe_tbl_reorg.mat stop_win_X_mfe_tbl_reorg;
save tAI_table_reorg.mat tAI_table_reorg;