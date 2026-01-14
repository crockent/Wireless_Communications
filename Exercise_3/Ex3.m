clear all; close all; clc;

%% Part A
K = 4;
N = 200;

%A1 - Channel creation (σωστό)
h = (randn(K+1,1) + 1j*randn(K+1,1))/sqrt(2);
energy = sum(abs(h).^2);
h_norm = h / sqrt(energy);  % Normalization
norm_energy = sum(abs(h_norm).^2);  % ~1

%A2
bits_all = randi([0 1], 2*N, 1);
symbol_package = bits_to_4qam(bits_all); 

%A3
n1 = 40;
n2 = 60;
training_symbols = symbol_package(n1:1:n2);  
L_train = length(training_symbols); % n2-n1+1 = 21

%A4
y = conv(symbol_package, h_norm);
%L_y = length(y);  % N+K

%A5       
Py = 2 * sum(abs(h_norm).^2); 

%A6 
SNR_db = 30;
SNR_linear = 10^(SNR_db/10);
N0 = Py / SNR_linear;
noise = sqrt(N0/2) * (randn(size(y)) + 1j*randn(size(y)));
y_n = y + noise;

%A7 
figure;
scatter(real(y), imag(y), 20, 'b', 'filled'); hold on;
scatter(real(y_n), imag(y_n), 20, 'r','filled');
grid on; axis equal;
xlabel('Re{y}'); ylabel('Im{y}');
legend('Without noise', 'With noise');
title(sprintf('Channel output, SNR = %d dB', SNR_db));
hold off;

%% Part B 
%B1-B2
y_train = y_n(n1+K : n2);  


col_input = symbol_package(n1+K : 1 : n2);  
row_input = symbol_package(n1+K : -1 : n1);  
X_train = toeplitz(col_input, row_input);  

%B3 - LS estimation
h_ls = (X_train' * X_train) \ (X_train' * y_train);
h_est_energy = sum(abs(h_ls).^2);

%B4 - Enhanced plot: Re/Im parts separately
figure;
plot_dim = 0:K;
subplot(2,1,1);
stem(plot_dim, real(h_norm), 'bo-', 'LineWidth', 1.2); hold on;
stem(plot_dim, real(h_ls), 'rx--', 'LineWidth', 1.2);
grid on; xlabel('Tap k'); ylabel('Re{h_k}');
legend('True', 'LS estimate'); title('Real part');
subplot(2,1,2);
stem(plot_dim, imag(h_norm), 'bo-', 'LineWidth', 1.2); hold on;
stem(plot_dim, imag(h_ls), 'ro--', 'LineWidth', 1.2);
grid on; xlabel('Tap k'); ylabel('Im{h_k}');
legend('True', 'LS estimate'); title('Imaginary part');


%% PART C 
%C1
delta = K;
f_zf = compute_zf_equalizer(h_ls, K, delta);

%C2
g = conv(h_norm,f_zf);

figure;

subplot(2,1,1);
stem(0:length(g)-1, real(g), 'bo-','LineWidth',1.2);
grid on;
xlabel('n');
ylabel('Re\{g^{ZF}\}');
title('Real part of \{h * f_{ZF}\}');

subplot(2,1,2);
stem(0:length(g)-1, imag(g), 'ro-','LineWidth',1.2);
grid on;
xlabel('n');
ylabel('Im\{g^{ZF}\}');
title('Imaginary part of \{h * f_{ZF}\}');


%C3
r_zf = conv(y_n,f_zf);


%C4
start_idx = delta+1 ;
end_ind = start_idx + N ;
s_n = r_zf(start_idx : end_ind);

%c5
constellation = [ 1+1i; 1-1i; -1+1i; -1-1i ];
figure;
scatter(real(s_n), imag(s_n), 20, 'b', 'filled'); hold on;
scatter(real(constellation), imag(constellation), 20, 'r', 'filled'); hold on;

% Add x=0 (vertical) and y=0 (horizontal) lines
plot([0 0], ylim, 'k--', 'LineWidth', 1); hold on;
plot(xlim, [0 0], 'k--', 'LineWidth', 1); hold on;

grid on; axis equal;
xlabel('Re{y}'); ylabel('Im{y}');
title('Properly truncate the equalizer output');
hold off;


% C6 - LS Equalizer
K_eq = length(f_zf) - 2; 
f_ls = compute_ls_equalizer(y_n, n1, n2, K_eq, delta,training_symbols);

% Effective response
gLS = conv(h_norm, f_ls);
figure;
subplot(2,1,1);
stem(0:length(gLS)-1, real(gLS), 'bo-', 'LineWidth', 1.2); grid on;
xlabel('n'); ylabel('Re(g_{LS})'); title('Real part of h * f_{LS}');
subplot(2,1,2);
stem(0:length(gLS)-1, imag(gLS), 'ro-', 'LineWidth', 1.2); grid on;
xlabel('n'); ylabel('Im(g_{LS})'); title('Imaginary part of h * f_{LS}');

% Equalized symbols (truncate to match transmitted symbols)
rLS = conv(y_n, f_ls); 
startidx_LS = delta + 1;
endind_LS = delta + N;
sn_LS = rLS(startidx_LS : endind_LS);
figure;
scatter(real(sn_LS), imag(sn_LS), 20, 'b', 'filled');
hold on;
scatter(real(constellation), imag(constellation), 50, 'r', 'filled');
plot([0 0], ylim, 'k--', 'LineWidth', 1);
plot(xlim, [0 0], 'k--', 'LineWidth', 1);
grid on; axis equal;
xlabel('Re(y)'); ylabel('Im(y)'); title('LS Equalizer Constellation');
hold off;


%C7
zf_error = zeros(4*K,1);
MMSE_error = zeros(4*K,1);
for delta = 0 : 1 : 4*K
    % Zero forcing eq
    f_zf = compute_zf_equalizer(h_ls, K, delta);
    r_zf = conv(y_n,f_zf);
    start_idx = delta+1 ;
    end_ind = start_idx + N ;
    s_n_zf = r_zf(start_idx : end_ind);
    
    zf_error(delta+1) = mean(abs(s_n_zf(n1:n2) - training_symbols(1:1:end)).^2);

    % MMSE eq
    K_eq = length(f_zf) - 2; 
    f_ls = compute_ls_equalizer(y_n, n1, n2, K_eq, delta,training_symbols);
    rLS = conv(y_n, f_ls); 
    startidx_LS = delta + 1;
    endind_LS = startidx_LS + N - 1;
    sn_LS = rLS(startidx_LS : endind_LS);
    
    MMSE_error(delta+1) =  mean(abs(sn_LS(n1:n2) - training_symbols(1:1:end)).^2);

end

figure;
plot(0 : 1 : 4*K,zf_error);
hold on;
plot(0 : 1 : 4*K,MMSE_error);
grid on;
hold off;
xlabel('\Delta');
ylabel('Eq. error');
legend('ZF equalizer','LS equalizer');
title('Equalization error');


%% D
clear all; close all; clc;
% D1
K = 4;
N = 200;
h = (randn(K+1,1) + 1j*randn(K+1,1))/sqrt(2);
energy = sum(abs(h).^2);
h_norm = h / sqrt(energy);  % Normalization
norm_energy = sum(abs(h_norm).^2);  % ~1

n1 = 40;
n2 = 60;

%D2
SNR_db = [2:2:30];
delta = K;
M_Pack = 300;  % Monte Carlo packets per SNR
BER_ZF = zeros(length(SNR_db), 1);
BER_LS = zeros(length(SNR_db), 1);
for i_SNR = 1:length(SNR_db)
    SNR_lin = 10^(SNR_db(i_SNR)/10);
    errors_ZF = 0; errors_LS = 0;
    L_train = n2 - n1 + 1;
    num_train_bits = 2 * L_train;
    
    for m = 1:M_Pack
        % Packet generation (unchanged)
        bits_all = randi([0 1], 2*N, 1);
        s_package = bits_to_4qam(bits_all);
        y_clean = conv(s_package, h_norm);
        Py = 2 * mean(abs(s_package).^2);
        N0 = Py / SNR_lin;
        noise = sqrt(N0/2) * (randn(size(y_clean)) + 1j*randn(size(y_clean)));
        y_n = y_clean + noise;
        
        % Training reference bits (fix 1: define correctly)
        bits_train_zf = bits_all((n1-1)*2+1 : n2*2);
        
        % Channel estimation (unchanged)
        training_symbols = s_package(n1:n2);
        y_train = y_n(n1+K : n2);
        col_input = s_package(n1+K : n2);
        row_input = s_package(n1+K : -1 : n1);
        X_train = toeplitz(col_input, row_input);
        h_ls = (X_train' * X_train) \ (X_train' * y_train);
        
        % ZF (fix 2: correct error count)
        f_zf = compute_zf_equalizer(h_ls, K, delta);
        r_zf = conv(y_n, f_zf);
        s_hat_zf = r_zf(delta+1 : delta + N);
        bits_hat_zf = qam4_to_bits(s_hat_zf);
        train_bits_hat_zf = bits_hat_zf((n1-1)*2+1 : n2*2);
        errors_ZF = errors_ZF + sum(train_bits_hat_zf ~= bits_train_zf);
        
        % LS (same fixes)
        K_eq = length(f_zf)-2;
        f_ls = compute_ls_equalizer(y_n, n1, n2, K_eq, delta, training_symbols);
        r_ls = conv(y_n, f_ls);
        s_hat_ls = r_ls(delta+1 : delta + N);
        bits_hat_ls = qam4_to_bits(s_hat_ls);
        train_bits_hat_ls = bits_hat_ls((n1-1)*2+1 : n2*2);
        errors_LS = errors_LS + sum(train_bits_hat_ls ~= bits_train_zf);
    end
    % Fix 3: normalize by total bits
    BER_ZF(i_SNR) = errors_ZF / (M_Pack * num_train_bits);
    BER_LS(i_SNR) = errors_LS / (M_Pack * num_train_bits);
end

%D4
h = sqrt(1/(2*(K+1))) * (randn(K+1,1) + 1i*randn(K+1,1));
energy = sum(abs(h).^2);
h_norm = h / sqrt(energy);  % Normalization
norm_energy = sum(abs(h_norm).^2);  % ~1
BER_ZF_new = zeros(length(SNR_db), 1);
BER_LS_new = zeros(length(SNR_db), 1);

for i_SNR = 1:length(SNR_db)
    SNR_lin = 10^(SNR_db(i_SNR)/10);
    errors_ZF = 0; errors_LS = 0;
    L_train = n2 - n1 + 1;
    num_train_bits = 2 * L_train;
    
    for m = 1:M_Pack
        % Packet generation (unchanged)
        bits_all = randi([0 1], 2*N, 1);
        s_package = bits_to_4qam(bits_all);
        y_clean = conv(s_package, h_norm);
        Py = 2 * mean(abs(s_package).^2);
        N0 = Py / SNR_lin;
        noise = sqrt(N0/2) * (randn(size(y_clean)) + 1j*randn(size(y_clean)));
        y_n = y_clean + noise;
        
        % Training reference bits (fix 1: define correctly)
        bits_train_zf = bits_all((n1-1)*2+1 : n2*2);
        
        % Channel estimation (unchanged)
        training_symbols = s_package(n1:n2);
        y_train = y_n(n1+K : n2);
        col_input = s_package(n1+K : n2);
        row_input = s_package(n1+K : -1 : n1);
        X_train = toeplitz(col_input, row_input);
        h_ls = (X_train' * X_train) \ (X_train' * y_train);
        
        % ZF (fix 2: correct error count)
        f_zf = compute_zf_equalizer(h_ls, K, delta);
        r_zf = conv(y_n, f_zf);
        s_hat_zf = r_zf(delta+1 : delta + N);
        bits_hat_zf = qam4_to_bits(s_hat_zf);
        train_bits_hat_zf = bits_hat_zf((n1-1)*2+1 : n2*2);
        errors_ZF = errors_ZF + sum(train_bits_hat_zf ~= bits_train_zf);
        
        % LS (same fixes)
        K_eq = length(f_zf)-2;
        f_ls = compute_ls_equalizer(y_n, n1, n2, K_eq, delta, training_symbols);
        r_ls = conv(y_n, f_ls);
        s_hat_ls = r_ls(delta+1 : delta + N);
        bits_hat_ls = qam4_to_bits(s_hat_ls);
        train_bits_hat_ls = bits_hat_ls((n1-1)*2+1 : n2*2);
        errors_LS = errors_LS + sum(train_bits_hat_ls ~= bits_train_zf);
    end
    % Fix 3: normalize by total bits
    BER_ZF_new(i_SNR) = errors_ZF / (M_Pack * num_train_bits);
    BER_LS_new(i_SNR) = errors_LS / (M_Pack * num_train_bits);
end

figure;
semilogy(SNR_db, BER_ZF, 'b-o', 'LineWidth', 2); 
hold on;
semilogy(SNR_db, BER_LS, 'r-s', 'LineWidth', 2);
semilogy(SNR_db, BER_ZF_new, 'r-s', 'LineWidth', 2);
semilogy(SNR_db, BER_LS_new, 'r-s', 'LineWidth', 2);
xlabel('SNR (dB)'); 
ylabel('Bit Error Probability'); 
legend('ZF', 'LS',"ZF_rand", "LS_rand"); 
grid on;
ylim([1e-5 1]);  % Optional: focus range
hold off;

