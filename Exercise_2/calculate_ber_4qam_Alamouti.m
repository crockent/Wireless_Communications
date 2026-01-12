function BER_Alamouti = calculate_ber_4qam_Alamouti(N, K, SNR_dB, h1_all, h2_all, bits_all, symbols_all)
    BER_Alamouti = zeros(size(SNR_dB));
    signal_power = 2;  % E{|s_k|^2} = 2 for 4-QAM
    
    for idx = 1:length(SNR_dB)
        total_bit_errors = 0;
        total_bits = 0;
        SNR_lin = 10^(SNR_dB(idx)/10);
        N0 = signal_power / SNR_lin;
        
        for k = 1:K
            % Use pre-generated channels, bits, and symbols
            h1 = h1_all(k);
            h2 = h2_all(k);
            bit_seq = bits_all(:, k);
            s_k = symbols_all(:, k);  % Full length N symbols
            
            % Alamouti block processing: split into odd/even symbols
            L = N/2;  % Number of Alamouti blocks (must be even N)
            s1 = s_k(1:2:end) / sqrt(2);  % First symbol each block, power normalized
            s2 = s_k(2:2:end) / sqrt(2);  % Second symbol each block, power normalized
            
            % Generate noise for L time slots at each receive antenna
            n1 = sqrt(N0/2) * (randn(L,1) + 1j*randn(L,1));
            n2 = sqrt(N0/2) * (randn(L,1) + 1j*randn(L,1));
            
            % Alamouti encoding + received signals (time slots 1 and 2)
            s1_star = conj(s1);
            s2_star = conj(s2);
            
            % Time slot 1: y1 = h1*s1 + h2*s2 + n1
            y1 = h1 * s1 + h2 * s2 + n1;
            % Time slot 2: y2 = h1*(-s2*) + h2*s1* + n2  
            y2 = h1 * (-s2_star) + h2 * s1_star + n2;
            
            % Alamouti combining (matched filter)
            h_norm = abs(h1)^2 + abs(h2)^2;
            R1 = (conj(h1) .* y1 + h2 .* conj(y2)) / sqrt(h_norm);
            R2 = (conj(h2) .* y1 - h1 .* conj(y2)) / sqrt(h_norm);
            
            % 4-QAM decision: map to ±1 ±1j constellation
            s1_hat = sign(real(R1)) + 1j * sign(imag(R1));
            s2_hat = sign(real(R2)) + 1j * sign(imag(R2));
            
            % Reconstruct full symbol sequence (length N)
            s_hat = zeros(N, 1);
            s_hat(1:2:end) = s1_hat;
            s_hat(2:2:end) = s2_hat;
            
            % Demap symbols to bits (consistent with other functions)
            num_bits = 2*N;
            bit_seq_hat = zeros(num_bits, 1);
            bit_seq_hat(1:2:end) = (real(s_hat) > 0);  % Real bit
            bit_seq_hat(2:2:end) = (imag(s_hat) > 0);  % Imag bit
            
            % Count bit errors
            bit_errors = sum(bit_seq ~= bit_seq_hat);
            total_bit_errors = total_bit_errors + bit_errors;
            total_bits = total_bits + num_bits;
        end
        
        BER_Alamouti(idx) = total_bit_errors / total_bits;
    end
end
