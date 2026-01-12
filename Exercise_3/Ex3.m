clear all,
close all,
clc;

%% Part A
K = 4;
N=100;

%A1

%Creation of complex channel
h = (randn(K+1,1) + 1j*randn(K+1,1))/sqrt(2);
%energy of complex channel
energy = sum(abs(h).^2);
%normalization
h_norm = h / sqrt(energy);
norm_energy = sum(abs(h_norm).^2);

%A2
bits_all = randi([0 1], 2*N,1);
symbol_package = zeros(N,1);
for k = 1:K
    symbol_package = bits_to_4qam(bits_all);
end


%A3
n1 = 40;
n2 = 60;

training_symbols = symbol_package(n1:n2);

%A4
y = conv(symbol_package,h_norm);

figure
grid on;
stem(y);
grid off;
legend("Output without noise")

%A6
P_y = 2;
SNR_db = 30;
SNR_linear = 10^(SNR_db/10);  
N0 = P_y / SNR_linear;

n=sqrt(N0/2) * (randn(size(y)) + 1j*randn(size(y)));

y_n = y+n;

%A7
figure;
grid on;
stem(y);
hold on;
stem(y_n);
grid off;
legend("Output without noise", "Output with noise")
hold off;

%% Part B
%B1
L_train = length(training_symbols);  % 21
P = L_train - K;  % 17
c = [training_symbols; zeros(K,1)];  
r = [training_symbols(1) zeros(1,K)]; 
X = toeplitz(c, r);                   %Toeplitz

%B3 LS estimation
X_train = X(1:L_train, :);            % 21 x 5 για LS
y_train = y_n(n1:n2);                    % Received training: 21 x 1
h_ls = (X_train' * X_train) \ (X_train' * y_train);  % LS estimate
h_est_energy = sum(abs(h_ls).^2);        % Έλεγχος ενέργειας ~1

%B4
figure;
stem(h_norm);
hold on;
stem(h_ls);
stem(h_norm - h_ls)
legend("Actual Channel", "Estimate Channel","Channel Dif")
hold off;

%% Part C
%C1
f_zf =  (X'*X)*X.'
