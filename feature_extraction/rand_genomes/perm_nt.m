function perm_seq = perm_nt(seq)
     perm_seq(randperm(length(seq))) = seq(1:end);
end

