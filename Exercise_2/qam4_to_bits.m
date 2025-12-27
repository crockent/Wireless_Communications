function bits = qam4_to_bits(symbols, const, bits_const)
    bits = zeros(2*length(symbols),1);
    for k = 1:length(symbols)
        [~,idx] = min(abs(symbols(k)-const).^2);
        bits(2*k-1:2*k) = bits_const(idx,:).';
    end
end