clc;
clear all;
close all;

% parameters
N1 = 1000; % for slow fading 
N2 = 200;  % for fast fading
Rh0 = 1;
b1 = 0.99; % slow
b2 = 0.1;  % fast

%% 1) channel h[k] 
total_ignored_samples = round(1/(1-b1)); % samples to ignore for slow fading 

% channels
h1 = rayleigh_flat_fading(N1, b1, Rh0);
h2 = rayleigh_flat_fading(N2, b2, Rh0);

% h1 with the correct samples 
h1_correct_samples = h1(total_ignored_samples+1:end);
N1_used = length(h1_correct_samples); % number of samples eventually used


%% 2) input s[k]
num_bits_1 = 2 * N1_used;
num_bits_2 = 2 * N2;

bit_seq_1 = (sign(randn(num_bits_1, 1)) + 1)/2;
bit_seq_2 = (sign(randn(num_bits_2, 1)) + 1)/2;

% from bits to 4-QAM symbols
s_k_1 = bits_to_4QAM(bit_seq_1);
s_k_2 = bits_to_4QAM(bit_seq_2);

%% 3) output r[k]
SNR_dB = 20; % change if needed
Es = 2;
sigma2_n = Es / (10^(SNR_dB/10));

% generate complex white Gaussian noise n[k] ~ CN(0, sigma2_n)
n_1 = sqrt(sigma2_n/2) * (randn(N1_used, 1) + 1i*randn(N1_used, 1));
n_2 = sqrt(sigma2_n/2) * (randn(N2, 1) + 1i*randn(N2, 1));

% channel output: r[k] = h[k]*s[k] + n[k]
r_1 = h1_correct_samples .* s_k_1 + n_1;
r_2 = h2 .* s_k_2 + n_2;

%% plots for questions (1)-(3):

% Figure : Channel comparison
figure;

% slow fading channel
subplot(2,1,1);
plot(1:N1_used, real(h1_correct_samples), 'b-', 'LineWidth', 1);
hold on;
plot(1:N1_used, imag(h1_correct_samples), 'r-', 'LineWidth', 1);
title(sprintf('Slow Fading Channel - b = %.2f', b1));
xlabel('Time k');
ylabel('Amplitude');
legend('Real', 'Imaginary', 'Location', 'best');
grid on;

% fast fading channel
subplot(2,1,2);
plot(1:N2, real(h2), 'b-', 'LineWidth', 1);
hold on;
plot(1:N2, imag(h2), 'r-', 'LineWidth', 1);
title(sprintf('Fast Fading Channel - b = %.2f', b2));
xlabel('Time k');
ylabel('Amplitude');
legend('Real', 'Imaginary', 'Location', 'best');
grid on;

% Figure : Channel magnitude and phase 
figure;

% magnitude - slow fading
subplot(2,2,1);
plot(1:N1_used, abs(h1_correct_samples), 'g-', 'LineWidth', 1);
title(sprintf('Slow Fading Magnitude - b = %.2f', b1));
xlabel('Time k');
ylabel('|h[k]|');
grid on;

% magnitude - fast fading
subplot(2,2,2);
plot(1:N2, abs(h2), 'm-', 'LineWidth', 1);
title(sprintf('Fast Fading Magnitude - b = %.2f', b2));
xlabel('Time k');
ylabel('|h[k]|');
grid on;

% phase - slow fading
subplot(2,2,3);
plot(1:N1_used, angle(h1_correct_samples), 'g-', 'LineWidth', 1);
title(sprintf('Slow Fading Phase - b = %.2f', b1));
xlabel('Time k');
ylabel('Phase (rad)');
grid on;

% phase - fast fading
subplot(2,2,4);
plot(1:N2, angle(h2), 'm-', 'LineWidth', 1);
title(sprintf('Fast Fading Phase - b = %.2f', b2));
xlabel('Time k');
ylabel('Phase (rad)');
grid on;

% Figure : 4-QAM constellations
figure;

subplot(1,2,1);
plot(real(s_k_1), imag(s_k_1), 'bo', 'MarkerSize', 5);
title('4-QAM Input (Slow Case)');
xlabel('In-Phase');
ylabel('Quadrature');
axis equal; grid on;
xlim([-2 2]); ylim([-2 2]);

subplot(1,2,2);
plot(real(s_k_2), imag(s_k_2), 'ro', 'MarkerSize', 5);
title('4-QAM Input (Fast Case)');
xlabel('In-Phase');
ylabel('Quadrature');
axis equal; grid on;
xlim([-2 2]); ylim([-2 2]);

% Figure : Output signals
figure;

subplot(2,1,1);
plot(1:N1_used, real(r_1), 'b-', 'LineWidth', 1);
hold on;
plot(1:N1_used, imag(r_1), 'r-', 'LineWidth', 1);
title(sprintf('Slow Fading Output - b = 0.99 (SNR=%ddB)', SNR_dB));
xlabel('Time k');
ylabel('Amplitude');
legend('Real', 'Imaginary', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(1:N2, real(r_2), 'b-', 'LineWidth', 1);
hold on;
plot(1:N2, imag(r_2), 'r-', 'LineWidth', 1);
title(sprintf('Fast Fading Output - b = 0.10 (SNR=%ddB)', SNR_dB));
xlabel('Time k');
ylabel('Amplitude');
legend('Real', 'Imaginary', 'Location', 'best');
grid on;

%% 4) ML detection with channel knowledge 

% rotate received signals: R = (h*/|h|) * r = |h| * s + Z
R_1 = (conj(h1_correct_samples) ./ abs(h1_correct_samples)) .* r_1;
R_2 = (conj(h2) ./ abs(h2)) .* r_2;

% ML detection: independent sign detection on real and imaginary parts
s_hat_1 = sign(real(R_1)) + 1j * sign(imag(R_1));
s_hat_2 = sign(real(R_2)) + 1j * sign(imag(R_2));

% calculate symbol error rates
ser_1 = sum(s_hat_1 ~= s_k_1) / length(s_k_1);
ser_2 = sum(s_hat_2 ~= s_k_2) / length(s_k_2);

fprintf('ML Detection with Perfect Channel Knowledge:\n');
fprintf('Slow fading (b=%.2f) SER: %.4f\n', b1, ser_1);
fprintf('Fast fading (b=%.2f) SER: %.4f\n', b2, ser_2);

%% plot for question (4):

% Figure : ML Detection results
figure;

subplot(1,2,1);
plot(real(s_k_1), imag(s_k_1), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
hold on;
plot(real(s_hat_1), imag(s_hat_1), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
title(sprintf('Slow Fading - b = 0.99 : Transmitted vs ML Detected \nSER=%.4f', ser_1));
xlabel('In-Phase');
ylabel('Quadrature');
axis equal; grid on;
xlim([-2 2]); ylim([-2 2]);
legend('Transmitted', 'ML Detected');

subplot(1,2,2);
plot(real(s_k_2), imag(s_k_2), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
hold on;
plot(real(s_hat_2), imag(s_hat_2), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
title(sprintf('Fast Fading- b = 0.1 : Transmitted vs ML Detected \nSER=%.4f', ser_2));
xlabel('In-Phase');
ylabel('Quadrature');
axis equal; grid on;
xlim([-2 2]); ylim([-2 2]);
legend('Transmitted', 'ML Detected');

%% 5) decoding errors

% absolute symbol error
err_abs_1 = abs(s_hat_1 - s_k_1);     
err_abs_2 = abs(s_hat_2 - s_k_2);

% Figure: Decoding errors and channel amplitude

figure;
subplot(2,1,1);
semilogy(abs(h1_correct_samples),'LineWidth',1);
title('Slow fading - b = 0.99 : |h[k]|'); xlabel('k'); ylabel('|h[k]|'); grid on;
subplot(2,1,2);
stem(err_abs_1,'filled');              
title('Slow fading - b = 0.99 : |s_{hat}[k]-s[k]|'); xlabel('k'); ylabel('|error|'); grid on;

figure;
subplot(2,1,1);
semilogy(abs(h2),'LineWidth',1);
title('Fast fading - b = 0.1 : |h[k]|'); xlabel('k'); ylabel('|h[k]|'); grid on;
subplot(2,1,2);
stem(err_abs_2,'filled');
title('Fast fading - b = 0.1 : |s_{hat}[k]-s[k]|'); xlabel('k'); ylabel('|error|'); grid on;

%% 7) flat fading - BER calculation, Monte Carlo

% parameters
SNR_dB = 0:2:30;
K = 5000;  
N1 = 1000;
N2 = 1000;

% theoretical curves
Es = 2;                    % 4-QAM energy        
SNR_lin  = 10.^(SNR_dB/10); % SNR = Es/N0
%N0 = Es / SNR_lin;   

BER_theoretical = 0.5 .* (1 - sqrt(SNR_lin ./ (2 + SNR_lin)));
BER_high_SNR    = 1 ./ (2 .* SNR_lin);

% Monte-Carlos
BER_1 = zeros(size(SNR_dB));
BER_2 = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    total_bit_errors_1 = 0;  total_bits_1 = 0;
    total_bit_errors_2 = 0;  total_bits_2 = 0;

    for k = 1:K
      % 1) channel h[k] 
        total_ignored_samples = round(1/(1-b1)); % samples to ignore for slow fading 

        % channels
        h1 = rayleigh_flat_fading(N1, b1, Rh0);
        h2 = rayleigh_flat_fading(N2, b2, Rh0);

        % h1 with the correct samples 
        h1_correct_samples = h1(total_ignored_samples+1:end);
        N1_used = length(h1_correct_samples); % number of samples eventually used

      % 2) input s[k]
        num_bits_1 = 2 * N1_used;
        num_bits_2 = 2 * N2;

        bit_seq_1 = (sign(randn(num_bits_1, 1)) + 1)/2;
        bit_seq_2 = (sign(randn(num_bits_2, 1)) + 1)/2;

        % from bits to 4-QAM symbols
        s_k_1 = bits_to_4QAM(bit_seq_1);
        s_k_2 = bits_to_4QAM(bit_seq_2);
        
      % 3) output r[k]
        Es = 2;
        sigma2_n = Es / (10^(SNR_dB(i)/10));

        % generate complex white Gaussian noise n[k] ~ CN(0, sigma2_n)
        n_1 = sqrt(sigma2_n/2) * (randn(N1_used, 1) + 1i*randn(N1_used, 1));
        n_2 = sqrt(sigma2_n/2) * (randn(N2, 1) + 1i*randn(N2, 1));

        % channel output: r[k] = h[k]*s[k] + n[k]
        r_1 = h1_correct_samples .* s_k_1 + n_1;
        r_2 = h2 .* s_k_2 + n_2;
        
      % 4) ML detection with channel knowledge 
        % rotate received signals: R = (h*/|h|) * r = |h| * s + Z
        R_1 = (conj(h1_correct_samples) ./ abs(h1_correct_samples)) .* r_1;
        R_2 = (conj(h2) ./ abs(h2)) .* r_2;

        % ML detection: independent sign detection on real and imaginary parts
        s_hat_1 = sign(real(R_1)) + 1j * sign(imag(R_1));
        s_hat_2 = sign(real(R_2)) + 1j * sign(imag(R_2));

        % calculate symbol error rates
        ser_1 = sum(s_hat_1 ~= s_k_1) / length(s_k_1);
        ser_2 = sum(s_hat_2 ~= s_k_2) / length(s_k_2);

      % 5) convert symbols to bits again
        % +1 -> 0,   -1 -> 1
        bit_seq_hat_1 = zeros(num_bits_1,1);
        bit_seq_hat_1(1:2:end) = (real(s_hat_1) < 0); 
        bit_seq_hat_1(2:2:end) = (imag(s_hat_1) < 0);
        
        bit_seq_hat_2 = zeros(num_bits_2,1);
        bit_seq_hat_2(1:2:end) = (real(s_hat_2) < 0);
        bit_seq_hat_2(2:2:end) = (imag(s_hat_2) < 0);
        
      % 6) count bit errors
        bit_errors_1 = sum(bit_seq_1 ~= bit_seq_hat_1);
        total_bit_errors_1 = total_bit_errors_1 + bit_errors_1;
        total_bits_1 = total_bits_1 + num_bits_1;

        bit_errors_2 = sum(bit_seq_2 ~= bit_seq_hat_2);
        total_bit_errors_2 = total_bit_errors_2 + bit_errors_2;
        total_bits_2 = total_bits_2 + num_bits_2;
    end

    BER_1(i) = total_bit_errors_1 / total_bits_1;
    BER_2(i) = total_bit_errors_2 / total_bits_2;
end

figure;
semilogy(SNR_dB, BER_1, 'bo-','LineWidth',1.5,'MarkerSize',5); hold on;
semilogy(SNR_dB, BER_2, 'ms-','LineWidth',1.5,'MarkerSize',5);
semilogy(SNR_dB, BER_theoretical, 'r--','LineWidth',2);
semilogy(SNR_dB, BER_high_SNR,   'g:','LineWidth',2);
grid on; xlabel('SNR (dB)'); ylabel('BER');
title('BER vs SNR (4-QAM, Rayleigh fading, coherent ML)');
legend('Sim slow (b=0.99)','Sim fast (b=0.10)', ...
       'Theory (Rayleigh)','High-SNR approx','Location','southwest');

%% 8) AWGN - BER calculation, Monte Carlo

% parameters
SNR_dB = 0:2:14;
K = 5000;  
N1 = 1000;

% theoretical curves
Es = 2;                     % 4-QAM energy        
SNR_lin  = 10.^(SNR_dB/10); % SNR = Es/N0   

BER_AWGN_theoretical = qfunc( sqrt(SNR_lin) );

% Monte-Carlos
BER_AWGN_1 = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    total_bit_errors_1 = 0;  total_bits_1 = 0;

    for k = 1:K
      % 1) AWGN channel  
        h1 = 1;

      % 2) input s[k]
        num_bits_1 = 2 * N1;
        bit_seq_1 = (sign(randn(num_bits_1, 1)) + 1)/2;
       
        % from bits to 4-QAM symbols
        s_k_1 = bits_to_4QAM(bit_seq_1);
        
      % 3) output r[k]
        Es = 2;
        sigma2_n = Es / (10^(SNR_dB(i)/10));

        % generate complex white Gaussian noise n[k] ~ CN(0, sigma2_n)
        n_1 = sqrt(sigma2_n/2) * (randn(N1, 1) + 1i*randn(N1, 1));

        % channel output: r[k] = s[k] + n[k]
        r_1 = h1 * s_k_1 + n_1; % h1 = 1
        
      % 4) ML detection with channel knowledge 
        % no rotation
        R_1 = r_1;

        % ML detection: independent sign detection on real and imaginary parts
        s_hat_1 = sign(real(R_1)) + 1j * sign(imag(R_1));

        % calculate symbol error rates
        ser_1 = sum(s_hat_1 ~= s_k_1) / length(s_k_1);

      % 5) convert symbols to bits again
        % +1 -> 0,   -1 -> 1
        bit_seq_hat_1 = zeros(num_bits_1,1);
        bit_seq_hat_1(1:2:end) = (real(s_hat_1) < 0); 
        bit_seq_hat_1(2:2:end) = (imag(s_hat_1) < 0);
        
      % 6) count bit errors
        bit_errors_1 = sum(bit_seq_1 ~= bit_seq_hat_1);
        total_bit_errors_1 = total_bit_errors_1 + bit_errors_1;
        total_bits_1 = total_bits_1 + num_bits_1;
    end

    BER_AWGN_1(i) = total_bit_errors_1 / total_bits_1;
end

figure;
semilogy(SNR_dB, BER_AWGN_1, 'bo-','LineWidth',1.5,'MarkerSize',5); hold on;
semilogy(SNR_dB, BER_AWGN_theoretical, 'r--','LineWidth',2);
grid on; xlabel('SNR (dB)'); ylabel('BER');
title('BER vs SNR (4-QAM, AWGN, coherent ML)');
legend(sprintf('AWGN sim (N=%d)',N1), ...
       'AWGN theory Q(sqrt(SNR))','Location','southwest');

%% plot for AWGN - Rayleigh comparison

% make same SNR for both
SNR_dB_rayleigh = 0:2:14; 
rayleigh_indices = 1:length(SNR_dB_rayleigh);

% keep SNR_dB of AWGN since it's smaller
figure;
semilogy(SNR_dB, BER_AWGN_1, 'bo-','LineWidth',1.2,'MarkerSize',5); hold on;
semilogy(SNR_dB, BER_1(rayleigh_indices), 'rs-','LineWidth',1.2,'MarkerSize',5);
semilogy(SNR_dB, BER_2(rayleigh_indices), 'ms-','LineWidth',1.2,'MarkerSize',5);

grid on; xlabel('SNR (dB)'); ylabel('BER');
title('BER comparison: AWGN vs Rayleigh (4-QAM, coherent ML)');
legend('AWGN sim', ...
       'Rayleigh sim (b=0.99)','Rayleigh sim (b=0.10)', ...
       'Location','southwest');
   
%% 9) block flat fading channels

% parameters
SNR_dB = 0:2:30;
K_blocks = 10000;  
N_block = 100;     

% theoretical curves 
Es = 2;                     % 4-QAM energy        
SNR_lin  = 10.^(SNR_dB/10); % SNR = Es/N0
BER_theoretical = 0.5 .* (1 - sqrt(SNR_lin ./ (2 + SNR_lin)));

% Monte-Carlo for block fading
BER_block = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    total_bit_errors_block = 0;
    total_bits_block = 0;

    for k_block = 1:K_blocks
      % 1) channel for this block (constant for entire block): h ~ CN(0,1)
        h_block = (randn(1) + 1i*randn(1)) / sqrt(2); % CN(0,1)
        
      % 2) input s[k] 
        num_bits_block = 2 * N_block;
        bit_seq_block = (sign(randn(num_bits_block, 1)) + 1)/2;
        
        % from bits to 4-QAM symbols
        s_k_block = bits_to_4QAM(bit_seq_block);
        
      % 3) output r[k] for this block
        Es = 2;
        sigma2_n = Es / (10^(SNR_dB(i)/10));
        
        % generate complex white Gaussian noise n[k] ~ CN(0, sigma2_n)
        n_block = sqrt(sigma2_n/2) * (randn(N_block, 1) + 1i*randn(N_block, 1));
        
        % channel output: r[k] = h_block * s[k] + n[k]
        r_block = h_block * s_k_block + n_block;
        
      % 4) ML detection with channel knowledge
        % rotate received signals: R = (h*/|h|) * r = |h| * s + Z
        R_block = (conj(h_block) / abs(h_block)) * r_block;
        
        % ML detection: independent sign detection on real and imaginary parts
        s_hat_block = sign(real(R_block)) + 1j * sign(imag(R_block));
        
      % 5) convert symbols to bits again
        % +1 -> 0,   -1 -> 1
        bit_seq_hat_block = zeros(num_bits_block, 1);
        bit_seq_hat_block(1:2:end) = (real(s_hat_block) < 0); 
        bit_seq_hat_block(2:2:end) = (imag(s_hat_block) < 0);
        
      % 6) count bit errors
        bit_errors_block = sum(bit_seq_block ~= bit_seq_hat_block);
        total_bit_errors_block = total_bit_errors_block + bit_errors_block;
        total_bits_block = total_bits_block + num_bits_block;
    end

    BER_block(i) = total_bit_errors_block / total_bits_block;
end

% plot for comparison
figure;
semilogy(SNR_dB, BER_block, 'k^-','LineWidth',1.5,'MarkerSize',6); hold on;
semilogy(SNR_dB, BER_1, 'bo-','LineWidth',1.5,'MarkerSize',5);
semilogy(SNR_dB, BER_2, 'ms-','LineWidth',1.5,'MarkerSize',5);
semilogy(SNR_dB, BER_theoretical, 'r--','LineWidth',2);
grid on; 
xlabel('SNR (dB)'); 
ylabel('BER');
title('BER vs SNR: Block Fading vs Time-Varying Fading (4-QAM)');
legend('Block fading (K=10000, N=100)', ...
       'Time-varying slow (b=0.99)', ...
       'Time-varying fast (b=0.10)', ...
       'Theory (Rayleigh)','Location','southwest');

