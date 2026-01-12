function BER = calculate_ber_4qam_tb(N, K, SNR_dB, h1_all, h2_all, bits_all, symbols_all)
    BER = zeros(size(SNR_dB));  
    signal_power = 2;
    for idx = 1:length(SNR_dB)
        total_bit_errors = 0;
        total_bits = 0;
        SNR_linear = 10^(SNR_dB(idx)/10);  
        N0 = signal_power / SNR_linear;
        for j = 1:K
            h1 = h1_all(j);
            h2 = h2_all(j);
            bits = bits_all(:, j);
            symbols = symbols_all(:, j);
            % Generate noise
            n = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));
            h_norm = sqrt(abs(h1)^2 + abs(h2)^2);
            R = h_norm * symbols + n;
            % Decision based on sign of real and imag parts
            s_hat = (real(R) >= 0) + 1i*(imag(R) >= 0);
            s_hat(s_hat == 0) = -1; % force -1 for zero
            estBits = zeros(2*N,1);
            estBits(1:2:end) = (real(s_hat) > 0);
            estBits(2:2:end) = (imag(s_hat) > 0);
            % Count bit errors
            bit_errors = sum(bits ~= estBits);
            total_bit_errors = total_bit_errors + bit_errors;
            total_bits = total_bits + length(bits);
        end
        BER(idx) = total_bit_errors / total_bits;
    end
end