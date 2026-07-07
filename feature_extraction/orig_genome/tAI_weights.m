% Sij
Sij = [ 0.9999 0.28 1 0 ; 1 1 0 1 ; 1 0 1 0.41 ; 0 1 0.68 1 ];

% tCGN map
% (http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/)
tCGN = containers.Map;
tCGN('AGC') = 11;
tCGN('TGC') = 5;
tCGN('GCC') = 18;
tCGN('CCC') = 2;
tCGN('TCC') = 3;
tCGN('AGG') = 2;
tCGN('TGG') = 10;
tCGN('AGT') = 11;
tCGN('CGT') = 1;
tCGN('TGT') = 4;
tCGN('AAC') = 14;
tCGN('CAC') = 2;
tCGN('TAC') = 2;
tCGN('AGA') = 11;
tCGN('CGA') = 1;
tCGN('TGA') = 3;
tCGN('GCT') = 4;
tCGN('ACG') = 6;
tCGN('CCG') = 1;
tCGN('CCT') = 1;
tCGN('TCT') = 11;
tCGN('GAG') = 1;
tCGN('TAG') = 3;
tCGN('CAA') = 10;
tCGN('TAA') = 7;
tCGN('GAA') = 10;
tCGN('GTT') = 10;
tCGN('CTT') = 14;
tCGN('TTT') = 7;
tCGN('GTC') = 16;
tCGN('CTC') = 2;
tCGN('TTC') = 14;
tCGN('GTG') = 7;
tCGN('CTG') = 1;
tCGN('TTG') = 9;
tCGN('AAT') = 13;
tCGN('TAT') = 2;
tCGN('CAT') = 5;
tCGN('GTA') = 8;
tCGN('GCA') = 4;
tCGN('CCA') = 6;

% tAI weights map
codons = [  "AAA" "AAC" "AAG" "AAT" ...
            "ACA" "ACC" "ACG" "ACT" ...
            "AGA" "AGC" "AGG" "AGT" ...
            "ATA" "ATC" "ATG" "ATT" ...
            "CAA" "CAC" "CAG" "CAT" ...
            "CCA" "CCC" "CCG" "CCT" ...
            "CGA" "CGC" "CGG" "CGT" ...
            "CTA" "CTC" "CTG" "CTT" ...
            "GAA" "GAC" "GAG" "GAT" ...
            "GCA" "GCC" "GCG" "GCT" ...
            "GGA" "GGC" "GGG" "GGT" ...
            "GTA" "GTC" "GTG" "GTT" ...
            "TAA" "TAC" "TAG" "TAT" ...
            "TCA" "TCC" "TCG" "TCT" ...
            "TGA" "TGC" "TGG" "TGT" ...
            "TTA" "TTC" "TTG" "TTT" ];
weights = zeros(1,length(codons));
for i = 1:1:length(codons)
    codon = char(codons(i));
    comp = seqrcomplement(codon);
    switch (codon(3))
        case 'A'
            j = 1;
        case 'C'
            j = 2;
        case 'G'
            j = 3;
        case 'T' %(U)
            j = 4;
    end
    tAI_w = 0;
    if tCGN.isKey(['A' comp(2:3)])
        tAI_w = tAI_w + tCGN(['A' comp(2:3)])*(1-Sij(1,j));
    end
    if tCGN.isKey(['C' comp(2:3)])
        tAI_w = tAI_w + tCGN(['C' comp(2:3)])*(1-Sij(2,j));
    end
    if tCGN.isKey(['G' comp(2:3)])
        tAI_w = tAI_w + tCGN(['G' comp(2:3)])*(1-Sij(3,j));
    end
    if tCGN.isKey(['T' comp(2:3)])
        tAI_w = tAI_w + tCGN(['T' comp(2:3)])*(1-Sij(4,j));
    end
    weights(i) = tAI_w;
end
max_w = max(weights);
notzeros = [];
for i = 1:1:length(codons)
    codon = codons(i);
    if weights(i) ~= 0
        weights(i) = weights(i)/max_w;
        notzeros = [notzeros weights(i)];
    end
end
tAI_w_map = containers.Map;
for i = 1:1:length(codons)
    codon = codons(i);
    if weights(i) ~= 0
        tAI_w_map(codon) = weights(i);
    else
        tAI_w_map(codon) = geomean(notzeros);
    end
end

save tAI_w_map.mat tAI_w_map;