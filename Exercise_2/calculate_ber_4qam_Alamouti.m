function BER = calculate_ber_4qam_Alamouti(N, K, SNR_dB, h1_all, h2_all, bits_all, symbols_all)
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
            n1 = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));
            n2 = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));
            % Received signals
            r1 = h1 * symbols + n1;
            r2 = h2 * symbols + n2;
            % Combine signals with conjugate weights
            r_total = conj(h1)*r1 + conj(h2)*r2;
            norm_h = abs(h1)^2 + abs(h2)^2;
            R = r_total / sqrt(norm_h);
            % Decision based on sign of real and imag parts
            s_hat = sign(real(R)) + 1i * sign(imag(R));
            % Reconstruct estimated bits from estimated symbols
            bit_seq_hat = zeros(2*N,1);
            bit_seq_hat(1:2:end) = (real(s_hat) > 0);
            bit_seq_hat(2:2:end) = (imag(s_hat) > 0);
            % Count bit errors
            bit_errors = sum(bits ~= bit_seq_hat);
            total_bit_errors = total_bit_errors + bit_errors;
            total_bits = total_bits + length(bits);
        end
        BER(idx) = total_bit_errors / total_bits;
    end
end