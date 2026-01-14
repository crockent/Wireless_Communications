function bits = qam4_to_bits(symbols)


    real_bits = (real(symbols) < 0);  %  Re<0 →1, else 0
    imag_bits = (imag(symbols) < 0);  %  Im<0 →1, else 0
    
    L = length(symbols);
    bits = zeros(2*L, 1);
    bits(1:2:end) = real_bits(:);
    bits(2:2:end) = imag_bits(:);

    
end
