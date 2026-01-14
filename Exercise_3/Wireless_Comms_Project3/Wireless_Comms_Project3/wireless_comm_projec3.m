clc;
clear all;
close all;

%% A1) 
K = 6;  

h = (randn(K+1,1) + 1i*randn(K+1,1)) / sqrt(2);
E_h = sum(abs(h).^2);

% normalization
h = h / sqrt(E_h);
E_h_norm = sum(abs(h).^2)

%% A2)
N = 200;                
num_bits = 2*N;         

bit_seq = (sign(randn(num_bits,1)) + 1)/2;

s = bits_to_4qam(bit_seq);   

constellation = [ 1+1i; 1-1i; -1+1i; -1-1i ];

%% A3) 
n1 = 30;        
n2 = 70;          

%% A4) 
y = conv(s, h);     


%% A5) 
% output power
Py_ex = mean(abs(y).^2)    
Py_th = 2 * sum(abs(h).^2) 

%% A6) 
SNR_dB  = 30;                      
SNRlin = 10^(SNR_dB/10);

sigma2_w = Py_th / SNRlin; 

w = sqrt(sigma2_w/2) * (randn(N+K,1) + 1i*randn(N+K,1));

yN = y + w;

%% A7)
figure;
scatter(real(y),  imag(y), 20, 'g', 'filled'); hold on;
scatter(real(yN), imag(yN), 20, 'b');
grid on;
axis equal;
xlabel('Re\{y_n\}');
ylabel('Im\{y_n\}');
legend('y_n (without noise)', 'y^N_n (with noise)');
title(['Channel output in the complex plane, SNR = ' num2str(SNR_dB) ' dB']);

%% B1-B3) 
first_col = s(n1+K : n2);          
first_row = s(n1+K : -1 : n1);    

A = toeplitz(first_col, first_row);   

y_train = yN(n1+K : n2); 

h_LS = (A' * A) \ (A' * y_train);

%% B4)
n_h = 0:K; % for plotting

figure;

subplot(2,1,1);
stem(n_h, real(h),   'bo-', 'LineWidth', 1.2); hold on;
stem(n_h, real(h_LS),'rx--','LineWidth', 1.2);
grid on;
xlabel('Tap index k');
ylabel('Re\{h_k\}');
legend('True h_k','LS estimate h_{LS,k}');
title('Real part of channel taps');

subplot(2,1,2);
stem(n_h, imag(h),   'bo-', 'LineWidth', 1.2); hold on;
stem(n_h, imag(h_LS),'rx--','LineWidth', 1.2);
grid on;
xlabel('Tap index k');
ylabel('Im\{h_k\}');
legend('True h_k','LS estimate h_{LS,k}');
title('Imaginary part of channel taps');

%% C1)
Lf = 4*K + 1; % equalizer length (k = 0,...,4K)
Lg = Lf + K;  % length of h*f

Delta = K;          

h_eq = h_LS;    

first_col_H = [h_eq; zeros(Lf-1,1)];     
first_row_H = [h_eq(1); zeros(Lf-1,1)];  
H_eq = toeplitz(first_col_H, first_row_H);  

eDelta = zeros(Lg,1);
eDelta(Delta + 1) = 1;     

% Zero-forcing equalizer
f_ZF = (H_eq' * H_eq) \ (H_eq' * eDelta);   

%% C2)
% overall impulse response with ZF 
g_ZF = conv(h, f_ZF);
n_g = 0:(Lg-1); % length = Lg

figure;

subplot(2,1,1);
stem(n_g, real(g_ZF), 'bo-','LineWidth',1.2);
grid on;
xlabel('n');
ylabel('Re\{g^{ZF}_n\}');
title('Real part of overall impulse response (h * f_{ZF})');

subplot(2,1,2);
stem(n_g, imag(g_ZF), 'bo-','LineWidth',1.2);
grid on;
xlabel('n');
ylabel('Im\{g^{ZF}_n\}');
title('Imaginary part of overall impulse response (h * f_{ZF})');

%% C3-C4)
z_ZF = conv(yN, f_ZF);
s_hat = z_ZF(Delta+1 : Delta+N); % estimates of s_n

%% C5) 
figure;
scatter(real(s_hat), imag(s_hat), 20, 'b'); hold on;
scatter(real(constellation), imag(constellation), 80, 'rx', 'LineWidth', 1.5);
grid on;
axis equal;
xlabel('Re(s\_hat)');
ylabel('Im(s\_hat)');
legend('s\_hat', '4-QAM constellation');
title(['Equalized symbols (ZF equalizer, SNR = ' num2str(SNR_dB) ' dB)']);

%% C6) LS
first_col_Y = yN(n1+Delta : n2+Delta);                     
first_row_Y = yN(n1+Delta : -1 : n1+Delta-(Lf-1));         
Y = toeplitz(first_col_Y, first_row_Y);                    

a_train = s(n1:n2);                                       

% LS equalizer
f_LS = (Y' * Y) \ (Y' * a_train);

% overall impulse response with LS 
g_LS = conv(h, f_LS);
n_g_LS = 0:(Lg-1);

figure;

subplot(2,1,1);
stem(n_g_LS, real(g_LS), 'bo-','LineWidth',1.2);
grid on;
xlabel('n');
ylabel('Re\{g^{LS}_n\}');
title('Real part of overall impulse response (h * f_{LS})');

subplot(2,1,2);
stem(n_g_LS, imag(g_LS), 'bo-','LineWidth',1.2);
grid on;
xlabel('n');
ylabel('Im\{g^{LS}_n\}');
title('Imaginary part of overall impulse response (h * f_{LS})');

% pass yN through LS equalizer and crop 
z_LS = conv(yN, f_LS);
s_hat_LS = z_LS(Delta+1 : Delta+N);

% plot equalized symbols for LS 
figure;
scatter(real(s_hat_LS), imag(s_hat_LS), 20, 'b'); hold on;
scatter(real(constellation), imag(constellation), 80, 'rx', 'LineWidth', 1.5);
grid on;
axis equal;
xlabel('Re(s\_hat\_LS)');
ylabel('Im(s\_hat\_LS)');
legend('s\_hat\_LS', '4-QAM constellation');
title(['Equalized symbols (LS equalizer, SNR = ' num2str(SNR_dB) ' dB)']);

%% C.7) 
Delta_values = 0:(4*K);
num_deltas   = length(Delta_values);

% H from h_LS (like C1)
first_col_H = [h_LS; zeros(Lf-1,1)];
first_row_H = [h_LS(1); zeros(Lf-1,1)];
H = toeplitz(first_col_H, first_row_H);  

% equalization error
eq_error_ZF = zeros(num_deltas,1);
eq_error_LS = zeros(num_deltas,1);

for idx_delta = 1:num_deltas
    Delta_current = Delta_values(idx_delta);

    % ===== ZF equalizer for this Delta =====
    eDelta = zeros(Lg,1);
    eDelta(Delta_current + 1) = 1;

    f_ZF_delta = (H' * H) \ (H' * eDelta);

    % estimates of s_n
    z_ZF_delta = conv(yN, f_ZF_delta);
    s_hat_ZF_delta = z_ZF_delta(Delta_current+1 : Delta_current+N);

    % eq. error on n = n1,...,n2
    eq_error_ZF(idx_delta) = mean(abs(s_hat_ZF_delta(n1:n2) - s(n1:n2)).^2);
 
    % ===== LS equalizer for this Delta =====
    first_col_Y = yN(n1 + Delta_current : n2 + Delta_current);
    first_row_Y = yN(n1 + Delta_current : -1 : n1 + Delta_current - (Lf-1));
    Y_delta = toeplitz(first_col_Y, first_row_Y);

    a_train = s(n1:n2);

    f_LS_delta = (Y_delta' * Y_delta) \ (Y_delta' * a_train);

    % estimates of s_n
    z_LS_delta = conv(yN, f_LS_delta);
    s_hat_LS_delta = z_LS_delta(Delta_current+1 : Delta_current+N);

    % eq. error on n = n1,...,n2
    eq_error_LS(idx_delta) = mean(abs(s_hat_LS_delta(n1:n2) - s(n1:n2)).^2);
end

% Plot 
figure;
plot(Delta_values, eq_error_ZF, 'bo-','LineWidth',1.2); hold on;
plot(Delta_values, eq_error_LS, 'ro-','LineWidth',1.2);
grid on;
xlabel('\Delta');
ylabel('Eq. error');
legend('ZF equalizer','LS equalizer');
title('Equalization error vs decision delay \Delta');

% best Delta 
[~, idxBestZF] = min(eq_error_ZF);
[~, idxBestLS] = min(eq_error_LS);

bestDelta_ZF = Delta_values(idxBestZF)
bestDelta_LS = Delta_values(idxBestLS)

% plots for g and s_hat for the best Delta (like C2 and C5/C6) 
n_g = 0:(Lg-1);

% ZF
eBest = zeros(Lg,1);
eBest(bestDelta_ZF + 1) = 1;

f_ZF_best = (H' * H) \ (H' * eBest);
g_ZF_best = conv(h, f_ZF_best);

z_ZF_best = conv(yN, f_ZF_best);
s_hat_ZF_best = z_ZF_best(bestDelta_ZF+1 : bestDelta_ZF+N);

% LS
first_col_Yb = yN(n1 + bestDelta_LS : n2 + bestDelta_LS);
first_row_Yb = yN(n1 + bestDelta_LS : -1 : n1 + bestDelta_LS - (Lf-1));
Y_best = toeplitz(first_col_Yb, first_row_Yb);

a_train = s(n1:n2);
f_LS_best = (Y_best' * Y_best) \ (Y_best' * a_train);
g_LS_best = conv(h, f_LS_best);

z_LS_best = conv(yN, f_LS_best);
s_hat_LS_best = z_LS_best(bestDelta_LS+1 : bestDelta_LS+N);

% figure 1
figure;
subplot(2,1,1);
stem(n_g, real(g_ZF_best), 'bo-','LineWidth',1.2); hold on;
stem(n_g, real(g_LS_best), 'gx--','LineWidth',1.2);
grid on;
xlabel('n'); ylabel('Re\{g_n\}');
legend(['ZF best \Delta=' num2str(bestDelta_ZF)], ['LS best \Delta=' num2str(bestDelta_LS)]);
title('Real part of overall impulse response (best \Delta)');

subplot(2,1,2);
stem(n_g, imag(g_ZF_best), 'bo-','LineWidth',1.2); hold on;
stem(n_g, imag(g_LS_best), 'gx--','LineWidth',1.2);
grid on;
xlabel('n'); ylabel('Im\{g_n\}');
legend(['ZF best \Delta=' num2str(bestDelta_ZF)], ['LS best \Delta=' num2str(bestDelta_LS)]);
title('Imaginary part of overall impulse response (best \Delta)');

%figure 2
figure;
scatter(real(s_hat_ZF_best), imag(s_hat_ZF_best), 20, 'b'); hold on;
scatter(real(s_hat_LS_best), imag(s_hat_LS_best), 20, 'g');
scatter(real(constellation), imag(constellation), 80, 'rx', 'LineWidth', 1.5);
grid on; axis equal;
xlabel('Re(s_hat)'); ylabel('Im(s_hat)');
legend('ZF (best \Delta)','LS (best \Delta)','4-QAM constellation');
title(['Equalized symbols for best \Delta (ZF=' num2str(bestDelta_ZF) ', LS=' num2str(bestDelta_LS) ')']);

%% D) 
SNRdB_vec = 0:2:30;
numSNR = length(SNRdB_vec);

M_packets = 500;

Lf = 4*K + 1;
Lg = Lf + K;
Delta = K;

% same training as before
s_train = s(n1:n2);

% LS matrix A is the same from B1

% data indices (exclude training)
idx_data = [1:n1-1, n2+1:N];
idx_bits_data = sort([2*idx_data-1, 2*idx_data]);

BER_ZF = zeros(numSNR,1);
BER_LS = zeros(numSNR,1);
BER_ZF_rand = zeros(numSNR,1);
BER_LS_rand = zeros(numSNR,1);

for iSNR = 1:numSNR
    SNR_dB = SNRdB_vec(iSNR);
    SNRlin = 10^(SNR_dB/10);

    for mode = 1:2
        % mode=1: fixed channel (D2)
        % mode=2: random channel per packet (D4)

        errZF = 0; errLS = 0; totBits = 0;

        for m = 1:M_packets

            % ===== choose channel for this packet =====
            if mode == 1
                h_use = h; % fixed, from A1 (already normalized)
            else
                h_use = sqrt(1/(2*(K+1))) * (randn(K+1,1) + 1i*randn(K+1,1)); %rand channel
            end

            Py = 2 * sum(abs(h_use).^2);     
            sigma2_w = Py / SNRlin;

            % build packet 
            bit_seq = (sign(randn(2*N,1)) + 1)/2;
            s_m = bits_to_4QAM(bit_seq);

            % insert training
            s_m(n1:n2) = s_train;

            y  = conv(s_m, h_use);           
            w  = sqrt(sigma2_w/2) * (randn(N+K,1) + 1i*randn(N+K,1));
            yN = y + w;

            % LS channel estimation from training 
            y_train = yN(n1+K : n2);        
            h_LS = (A' * A) \ (A' * y_train);

            % ===== ZF equalizer =====
            first_col_H = [h_LS; zeros(Lf-1,1)];
            first_row_H = [h_LS(1); zeros(Lf-1,1)];
            H_eq = toeplitz(first_col_H, first_row_H);   

            eDelta = zeros(Lg,1);
            eDelta(Delta+1) = 1;

            f_ZF = (H_eq' * H_eq) \ (H_eq' * eDelta);   

            % ===== LS equalizer =====
            y_col = yN(n1+Delta : n2+Delta);  
            y_row = yN(n1+Delta : -1 : n1+Delta-(Lf-1));
            Y = toeplitz(y_col, y_row);

            f_LS = (Y' * Y) \ (Y' * s_train);

            % equalize + crop 
            z_ZF = conv(yN, f_ZF);
            s_hat_ZF = z_ZF(Delta+1 : Delta+N);

            z_LS = conv(yN, f_LS);
            s_hat_LS = z_LS(Delta+1 : Delta+N);

            % hard decisions 
            bits_hat_ZF = zeros(2*N,1);
            bits_hat_ZF(1:2:end) = real(s_hat_ZF) < 0;
            bits_hat_ZF(2:2:end) = imag(s_hat_ZF) < 0;

            bits_hat_LS = zeros(2*N,1);
            bits_hat_LS(1:2:end) = real(s_hat_LS) < 0;
            bits_hat_LS(2:2:end) = imag(s_hat_LS) < 0;

            % count errors only on data bits 
            errZF = errZF + sum(bits_hat_ZF(idx_bits_data) ~= bit_seq(idx_bits_data));
            errLS = errLS + sum(bits_hat_LS(idx_bits_data) ~= bit_seq(idx_bits_data));
            totBits = totBits + length(idx_bits_data);
        end

        berZF = errZF / totBits;
        berLS = errLS / totBits;

        if mode == 1
            BER_ZF(iSNR) = berZF;
            BER_LS(iSNR) = berLS;
        else
            BER_ZF_rand(iSNR) = berZF;
            BER_LS_rand(iSNR) = berLS;
        end
    end

end

figure;
semilogy(SNRdB_vec, BER_ZF,      '-',  'LineWidth', 1.5); hold on;
semilogy(SNRdB_vec, BER_LS,      '-',  'LineWidth', 1.5);
semilogy(SNRdB_vec, BER_ZF_rand, '--', 'LineWidth', 1.5);
semilogy(SNRdB_vec, BER_LS_rand, '--', 'LineWidth', 1.5);
grid on;
xlabel('SNR in dB'); ylabel('Bit Error Rate');
legend('BER\_ZF','BER\_LS','BER\_ZF\_rand','BER\_LS\_rand','Location','northeast');
title('BER performance (ZF vs LS) - fixed and random channel');

