clear;
yeast_parameters_table = readtable("yeast_parameters_table.xls");
SGD_all_ORFs_5prime_UTRs = fastaread("SGD_all_ORFs_5prime_UTRs.fsa");
SGD_all_ORFs_3prime_UTRs = fastaread("SGD_all_ORFs_3prime_UTRs.fsa");
orf_coding_all = fastaread("orf_coding_all_R64-3-1_20210421.fasta");
Y_genes = cell(height(orf_coding_all),5);
gene_index = 0;

for i = 1:height(orf_coding_all)
    seq = orf_coding_all(i).Sequence;
    if mod(length(seq),3)==0
        gene_index = gene_index+1;
        Y_genes(gene_index,1) = textscan(orf_coding_all(i).Header,'%s*');
        Y_genes(gene_index,2) = {seq};
    end
end
for i = 1:height(SGD_all_ORFs_5prime_UTRs)
    name = textscan(SGD_all_ORFs_5prime_UTRs(i).Header,...
        'sacCer3_ct_PelechanoonlybasedUTRs_1122_%s*','delimiter','_');
    index = find(strcmp([Y_genes{:,1}],name{1}));
    if ~isempty(index)
        Y_genes(index,3) = {SGD_all_ORFs_5prime_UTRs(i).Sequence(1:end-1)};
    end
end
for i = 1:height(SGD_all_ORFs_3prime_UTRs)
    name = textscan(SGD_all_ORFs_3prime_UTRs(i).Header,...
        'sacCer3_ct_Pelechanoonlybased3primeUTRs_3950_%s*','delimiter','_');
    index = find(strcmp([Y_genes{:,1}],name{1}));
    if ~isempty(index)
        Y_genes(index,4) = {SGD_all_ORFs_3prime_UTRs(i).Sequence(2:end)};
    end
end
Y_genes = Y_genes(1:gene_index,:);

for i = 1:height(yeast_parameters_table)
    name = yeast_parameters_table.ORF(i);
    index = find(strcmp([Y_genes{:,1}],name{1}));
    if ~isnan(yeast_parameters_table.PA1(i))
        Y_genes(index,5) = {yeast_parameters_table.PA1(i)};
    end
end

save Y_genes.mat Y_genes;