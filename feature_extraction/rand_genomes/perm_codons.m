function perm_seq = perm_codons(seq,n)
    aminos = ['A' 'R' 'N' 'D' 'C' 'Q' 'E' 'G' 'H' 'I' 'L' 'K' 'M' 'F' 'P'...
        'S' 'T' 'W' 'Y' 'V' 'B' 'Z' 'X'];
    aa = nt2aa(seq);
    perm_seq = cell(n,1);
    for i = 1:length(aminos)
        aa_ind = strfind(aa,aminos(i));
        aa_perm = [];
        if ~isempty(aa_ind)
            for j = 1:n
                aa_perm(1:length(aa_ind)) = aa_ind(randperm(length(aa_ind)));
                perm_seq{j}(aa_perm.*3-2) = seq(aa_ind.*3-2);
                perm_seq{j}(aa_perm.*3-1) = seq(aa_ind.*3-1);
                perm_seq{j}(aa_perm.*3) = seq(aa_ind.*3);
            end
        end
    end
    
    %handle stop codons
    stop_ind = strfind(aa,'*');
    for j = 1:n
        perm_seq{j}(stop_ind.*3-2) = seq(stop_ind.*3-2);
        perm_seq{j}(stop_ind.*3-1) = seq(stop_ind.*3-1);
        perm_seq{j}(stop_ind.*3) = seq(stop_ind.*3);
    end
end

