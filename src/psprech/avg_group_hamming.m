function agh = avg_group_hamming(as,as0)
    agh = mean(sum(xor(as,as0)));
end
