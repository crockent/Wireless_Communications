close all;
clear all;
clc
N = 200;
b_values = [0.15, 0.30];
noise_power = 0.1; 
h_all = cell(length(b_values)+1, 1); % Cell array to store h for each b

for ib = 1:length(b_values)
    b = b_values(ib);
    sigma2 = 1 - b^2; % Ensure E{|h[k]|^2} = 1
    h = zeros(1, N);
    h(1) = 1;
    e = sqrt(sigma2/2) * (randn(1, N) + 1i*randn(1, N));
    for k = 2:N
        h(k) = b*h(k-1) + e(k);
    end
    
    h_all{ib} = h; % Save h before plotting
    figure; % new figure for each b
    
    x_axis = 1:N;
    subplot(2,1,1);
    hold on;
    plot(x_axis, real(h), 'DisplayName', 'Real part');
    plot(x_axis, imag(h), 'DisplayName', 'Imag part');
    title(['N=', num2str(N), ', b=', num2str(b), ' - Real & Imaginary']);
    xlabel('k');
    ylabel('Amplitude');
    legend show;
    hold off;
    
    subplot(2,1,2);
    plot(x_axis, abs(h), 'DisplayName', 'Magnitude', 'LineWidth', 1.5);
    title(['N=', num2str(N), ', r=', num2str(b), ' - Magnitude']);
    xlabel('k');
    ylabel('Amplitude');
    legend show;
end

% Special case: N=1000, b=0.99 
N = 1000;
b = 0.99;
sigma2 = 1 - b^2;
h = zeros(1, N);
h(1) = 1;
e = sqrt(sigma2/2) * (randn(1, N) + 1i*randn(1, N));
for k = 2:N
    h(k) = b*h(k-1) + e(k);
end
h_all{3} = h;
ignore_samples = ceil(1/(1-b));
h_plot = h(ignore_samples+1:end);
x_axis = (ignore_samples+1):N;

figure;
subplot(2,1,1);
hold on;
plot(x_axis, real(h_plot), 'DisplayName', 'Real part');
plot(x_axis, imag(h_plot), 'DisplayName', 'Imag part');
title(['N=', num2str(N), ', b=', num2str(b), ' - Real & Imaginary']);
xlabel('k');
ylabel('Amplitude');
legend show;
hold off;

subplot(2,1,2);
plot(x_axis, abs(h_plot), 'DisplayName', 'Magnitude', 'LineWidth', 1.5);
title(['N=', num2str(N), ', r=', num2str(b), ' - Magnitude']);
xlabel('k');
ylabel('Amplitude');
legend show;


%% 2- 3 
%--------- PARAMETERS ---------
SNR_dB = 20;          
Es = 2;               

%{
    r(1) is for N=200 b=0.15
    r(2) is for N=200 b=0.3
    r(3) is for N=1000 b=0.99
%}
r = cell(3, 1);
s = cell(3, 1);
noise_power =  Es / (10^(SNR_dB/10));
for i = 1:length(h_all)
    % Δημιουργία QAM συμβόλων
    s{i} = randsrc(1, length(h_all{i}), [1 -1]) + 1i * randsrc(1, length(h_all{i}), [1 -1]);
    % Λευκό Gaussian θόρυβο
    n = sqrt(noise_power/2) * (randn(1, length(h_all{i})) + 1i * randn(1, length(h_all{i})));
    % Έξοδος καναλιού: r[k] = h[k] * s[k] + n[k]
    r{i} = h_all{i} .* s{i} + n;
end

for i = 1:length(r)
    figure;
    plot(real(r{i}), 'b', 'DisplayName', 'Real part');
    hold on;
    plot(imag(r{i}), 'r', 'DisplayName', 'Imag part');
    hold off;
    grid on;
    if i <= 2
        title(['Received Signal r for N=200, b=', num2str(b_values(i))]);
    else
        title(['Received Signal r for N=1000, b=0.99']);
    end
    xlabel('Sample Index k');
    ylabel('Amplitude');
    legend show;
end

%% 5
qam_syms = [1+1j, 1-1j, -1+1j, -1-1j];
%------- ML Detection--------
ML_hat = cell(3,1); % Εκτιμήσεις για κάθε περίπτωση καναλιού
for i = 1:length(r)
    r_curr = r{i};
    h_curr = h_all{i};
    ML_detected = zeros(1, length(r_curr));
    for k = 1:length(r_curr)
        dists = abs(r_curr(k) - h_curr(k) * qam_syms);
        [~, idx] = min(dists);
        ML_detected(k) = qam_syms(idx);
    end
    ML_hat{i} = ML_detected;
end

%% 5 
threshold = 1/sqrt(SNR_dB);
for idx = 1:3
    s_true = s{idx};
    s_est = ML_hat{idx};
    errors = (s_true ~= s_est);    % 1 εκεί όπου υπάρχει σφάλμα
    h_amp = abs(h_all{idx});       % Πλάτος καναλιού |h[k]|
    
    figure;
    subplot(2,1,1);
    semilogy(h_amp, 'b', 'DisplayName', '|h[k]|'); % ημιλογαριθμική κλίμακα
    hold on;
    stem(find(errors), h_amp(errors), 'r', 'filled', 'DisplayName', 'Error positions');
    yline(threshold, 'k--', 'LineWidth', 1.5, 'DisplayName', 'threshold');
    xlabel('Sample Index');
    ylabel('|h[k]|');
    title(['Channel amplitude (semilogy) and error positions, channel ', num2str(idx)]);
    legend('show');
    hold off;
    
    subplot(2,1,2);
    stem(abs(errors), 'r', 'DisplayName', 'Decoding error (1=error)');
    xlabel('Sample Index');
    ylabel('Error');
    title(['Decoded symbol errors, channel ', num2str(idx)]);
    legend('show');
end


%%7
% Parameters
SNR_dB = 0:2:30;
SNR_lin = 10.^(SNR_dB/10);
K = 5000;           % Number of frames
N = 1000;           % Symbols per frame
b_values = [0.15 0.3 0.99];
Es = 2;             % QAM symbol energy

BER_sim = zeros(length(b_values), length(SNR_dB));

for b_idx = 1:length(b_values)
    b = b_values(b_idx);
    sigma2 = 1 - b^2; % So that E[|h|^2] = 1
    for snr_idx = 1:length(SNR_dB)
        snr_curr = SNR_dB(snr_idx);
        sigma2_n = Es / (10^(snr_curr/10));
        total_bit_errors = 0;
        total_bits = 0;

        for frame = 1:K
            % Channel generation (AR(1) Rayleigh)
            h = zeros(1, N);
            h(1) = 1;
            noise_h = sqrt(sigma2/2) * (randn(1, N) + 1i*randn(1, N));
            for n = 2:N
                h(n) = b * h(n-1) + noise_h(n);
            end
            % Optionally ignore transient samples for b ≈ 1
            if b > 0.9
                ignore_samples = round(1/(1-b));
                h_used = h(ignore_samples+1:end);
                n_sym = length(h_used);
            else
                h_used = h;
                n_sym = N;
            end

            % Symbol and noise
            bits = randi([0 1], 2*n_sym, 1); % 2 bits per symbol
            s = (2*bits(1:2:end)-1) + 1i*(2*bits(2:2:end)-1); % 4-QAM mapping
            n = sqrt(sigma2_n/2) * (randn(n_sym, 1) + 1i*randn(n_sym, 1));

            % Received signal
            r = h_used.' .* s + n;

            % Coherent detection
            r_eq = (conj(h_used)./abs(h_used)).' .* r;

            % ML decision
            s_hat = sign(real(r_eq)) + 1i*sign(imag(r_eq));

            % Map back to bits
            bits_hat = zeros(2*n_sym,1);
            bits_hat(1:2:end) = real(s_hat) > 0;
            bits_hat(2:2:end) = imag(s_hat) > 0;

            % Count errors
            total_bit_errors = total_bit_errors + sum(bits_hat ~= bits);
            total_bits = total_bits + 2*n_sym;
        end

        BER_sim(b_idx, snr_idx) = total_bit_errors / total_bits;
    end
end


BER_theoretical = 0.5 .* (1 - sqrt(SNR_lin ./ (2 + SNR_lin)));
BER_high_SNR = 1 ./ (2 .* SNR_lin);

figure;
hold on;
semilogy(SNR_dB, BER_sim(1,:), 'ro-',   'LineWidth',2, 'MarkerSize',8); % b = 0.15, red circles solid
semilogy(SNR_dB, BER_sim(2,:), 'bs--',  'LineWidth',2, 'MarkerSize',8); % b = 0.3, blue squares dashed
semilogy(SNR_dB, BER_sim(3,:), 'g^:',   'LineWidth',2, 'MarkerSize',8); % b = 0.99, green triangles dotted

semilogy(SNR_dB, BER_theoretical, 'r--','LineWidth',2);
semilogy(SNR_dB, BER_high_SNR,'p-','LineWidth',2);
grid on; xlabel('SNR (dB)'); ylabel('BER');
title('BER vs SNR');
yscale log
hold off;


%% 8
% Parameters
SNR_dB = 0:2:14;
SNR_lin = 10.^(SNR_dB/10);
K = 5000;
N = 1000; % symbols (same N as Rayleigh)
Es = 2;

% Pre-allocate results
BER_AWGN_sim = zeros(1, length(SNR_dB));

for snr_idx = 1:length(SNR_dB)
    sigma2_n = Es / SNR_lin(snr_idx);
    total_bit_errors = 0;
    total_bits = 0;

    for frame = 1:K
        % Symbol generation (same as Rayleigh)
        bits = randi([0 1], 2*N, 1); % 2 bits/symbol (4-QAM)
        s = (2*bits(1:2:end)-1) + 1i*(2*bits(2:2:end)-1); % 4-QAM mapping

        % AWGN channel (no fading): h=1
        n = sqrt(sigma2_n/2) * (randn(N, 1) + 1i*randn(N, 1));
        r = s + n;

        % ML detection: sign on real and imag
        s_hat = sign(real(r)) + 1i*sign(imag(r));
        bits_hat = zeros(2*N,1);
        bits_hat(1:2:end) = real(s_hat) > 0;
        bits_hat(2:2:end) = imag(s_hat) > 0;

        % Count errors
        total_bit_errors = total_bit_errors + sum(bits_hat ~= bits);
        total_bits = total_bits + 2*N;
    end

    BER_AWGN_sim(snr_idx) = total_bit_errors / total_bits;
end

% Theoretical curve for AWGN (4-QAM, Gray coding, see note)
BER_AWGN_theory = qfunc(sqrt(SNR_lin));

figure; hold on;
% AWGN sim/theory
semilogy(SNR_dB, BER_AWGN_sim,    'r-',  'LineWidth',2, 'MarkerSize',8);
semilogy(SNR_dB, BER_AWGN_theory, 'b--',  'LineWidth',2);
yscale log
grid on; xlabel('SNR (dB)'); ylabel('BER');
title('BER Comparison: AWGN vs Rayleigh (4-QAM)');
legend('AWGN Sim', 'AWGN Theory', ...
       'Location', 'southwest');
hold off;

%% 9
% --- PARAMETERS ---
K = 10000;          % Number of packets/blocks
N = 100;            % Number of symbols per packet/block
Es = 2;             % QAM symbol energy
SNR_dB = 0:2:30;    % SNR range
SNR_lin = 10.^(SNR_dB/10);

% --- THEORETICAL BER (flat Rayleigh) ---
BER_theoretical = 0.5 .* (1 - sqrt(SNR_lin ./ (2 + SNR_lin)));

% --- BLOCK FLAT FADING SIMULATION ---
BER_block_flat = zeros(size(SNR_dB));
for snr_idx = 1:length(SNR_dB)
    total_bit_errors = 0;
    total_bits = 0;
    sigma2_n = Es / SNR_lin(snr_idx);
    for pkt = 1:K
        % 1. Channel coefficient for whole block
        h_block = (randn(1) + 1i*randn(1)) / sqrt(2); % CN(0,1)
        % 2. Data bits + 4-QAM modulation
        bits = randi([0 1], 2*N, 1); % 2 bits/symbol
        s = (2*bits(1:2:end)-1) + 1i*(2*bits(2:2:end)-1); % 4-QAM
        % 3. Noise
        n = sqrt(sigma2_n/2) * (randn(N,1) + 1i*randn(N,1));
        % 4. Received signal
        r = h_block * s + n;
        % 5. Equalization (coherent receiver)
        r_eq = (conj(h_block)/abs(h_block)) * r;
        % 6. Decision
        s_hat = sign(real(r_eq)) + 1i * sign(imag(r_eq));
        bits_hat = zeros(2*N,1);
        bits_hat(1:2:end) = real(s_hat) > 0;
        bits_hat(2:2:end) = imag(s_hat) > 0;
        % 7. Count errors
        total_bit_errors = total_bit_errors + sum(bits_hat ~= bits);
        total_bits = total_bits + 2*N;
    end
    BER_block_flat(snr_idx) = total_bit_errors / total_bits;
end

% --- TIME-VARYING RAYLEIGH FADING (symbol-wise) ---
BER_timevarying = zeros(size(SNR_dB));
for snr_idx = 1:length(SNR_dB)
    total_bit_errors = 0;
    total_bits = 0;
    sigma2_n = Es / SNR_lin(snr_idx);
    for pkt = 1:K
        % 1. Channel coefficient for each symbol (independent)
        h_sym = (randn(N,1) + 1i*randn(N,1)) / sqrt(2);
        % 2. Data bits + 4-QAM modulation
        bits = randi([0 1], 2*N, 1);
        s = (2*bits(1:2:end)-1) + 1i*(2*bits(2:2:end)-1);
        % 3. Noise
        n = sqrt(sigma2_n/2) * (randn(N,1) + 1i*randn(N,1));
        % 4. Received signal
        r = h_sym .* s + n;
        % 5. Equalization (coherent receiver)
        r_eq = (conj(h_sym)./abs(h_sym)) .* r;
        % 6. Decision
        s_hat = sign(real(r_eq)) + 1i * sign(imag(r_eq));
        bits_hat = zeros(2*N,1);
        bits_hat(1:2:end) = real(s_hat) > 0;
        bits_hat(2:2:end) = imag(s_hat) > 0;
        % 7. Count errors
        total_bit_errors = total_bit_errors + sum(bits_hat ~= bits);
        total_bits = total_bits + 2*N;
    end
    BER_timevarying(snr_idx) = total_bit_errors / total_bits;
end

% --- PLOTTING ---
figure;
semilogy(SNR_dB, BER_block_flat, 'mx-','LineWidth',2,'MarkerSize',7); hold on;
semilogy(SNR_dB, BER_timevarying, 'bo-','LineWidth',2,'MarkerSize',7);
semilogy(SNR_dB, BER_theoretical, 'r--','LineWidth',2);
grid on; 
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR: Block Flat Fading vs Time-varying Rayleigh Fading');
legend('Block flat fading (K=10000,N=100)', ...
       'Time-varying Rayleigh', ...
       'Theoretical Rayleigh','Location','southwest');
yscale log;
hold off;





