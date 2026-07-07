clear;
load('Y_genes.mat');

chr.chrI = fastaread("chr01.fsa").Sequence; chr.chrII = fastaread("chr02.fsa").Sequence;
chr.chrIII = fastaread("chr03.fsa").Sequence; chr.chrIV = fastaread("chr04.fsa").Sequence;
chr.chrV = fastaread("chr05.fsa").Sequence; chr.chrVI = fastaread("chr06.fsa").Sequence;
chr.chrVII = fastaread("chr07.fsa").Sequence; chr.chrVIII = fastaread("chr08.fsa").Sequence;
chr.chrIX = fastaread("chr09.fsa").Sequence; chr.chrX = fastaread("chr10.fsa").Sequence;
chr.chrXI = fastaread("chr11.fsa").Sequence; chr.chrXII = fastaread("chr12.fsa").Sequence;
chr.chrXIII = fastaread("chr13.fsa").Sequence; chr.chrXIV = fastaread("chr14.fsa").Sequence;
chr.chrXV = fastaread("chr15.fsa").Sequence; chr.chrXVI = fastaread("chr16.fsa").Sequence;

fileID = fopen('3prime_UTR_summary_updated.txt');
orig_utrs = textscan(fileID,'%s %s %d %d %*s %*s','HeaderLines',1);

Y_genes_orig3utrs = Y_genes;
Y_genes_orig3utrs(:,3) = cell(height(Y_genes_orig3utrs),1);

for i = 1:height(orig_utrs{1})
    name = orig_utrs{1}{i};
    curr_chr = orig_utrs{2}{i};
    utr3_start = orig_utrs{3}(i);
    utr3_end = orig_utrs{4}(i);
    
    index = find(strcmp([Y_genes{:,1}],name));

    if contains(name,'W')
        Y_genes_orig3utrs{index,3} = chr.(curr_chr)(utr3_start:utr3_end);
    else % C
        Y_genes_orig3utrs{index,3} = seqrcomplement(chr.(curr_chr)(utr3_start-1:utr3_end-1));
    end
end

orf_coding_all = fastaread("orf_coding_all_R64-3-1_20210421.fasta");

for i = 1:height(orf_coding_all)
    name = textscan(orf_coding_all(i).Header,'%s*');
    index = find(strcmp([Y_genes_orig3utrs{:,1}],name{1}));
    if ~isempty(index)
        if isempty(Y_genes_orig3utrs{index,3}) % missing 3'UTR
            curr_chr = strcat('chr',string(textscan(orf_coding_all(i).Header,'%*s %*s %*s %*s %s*')));
            orf_end = textscan(orf_coding_all(i).Header,'%*s %*s %*s %*s %*s %*s %*d- %d*');
            orf_start = textscan(orf_coding_all(i).Header,'%*s %*s %*s %*s %*s %*s %d*- %*d');
            if ~strcmp(curr_chr,'chrMito')
                if contains(name{1},'W')
                    Y_genes_orig3utrs{index,3} = chr.(string(curr_chr))(orf_end{1}+1:orf_end{1}+173);
                else % C
                    Y_genes_orig3utrs{index,3} = seqrcomplement(chr.(string(curr_chr))...
                        (orf_end{1}-173:orf_end{1}-1));
                end
            end
        end
    end
end

save Y_genes_orig3utrs.mat Y_genes_orig3utrs;
