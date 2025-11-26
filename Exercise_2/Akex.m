clc;
clear all;
close all;

N = 200;

%% 1)
% CN(0,1)
h1 = (randn(1,1) + 1i*randn(1,1)) / sqrt(2);   
h2 = (randn(1,1) + 1i*randn(1,1)) / sqrt(2);   

% block fading
h1_block = h1*ones(N, 1);   
h2_block = h2*ones(N, 1);   

%% 2)
num_bits = 2*N;
bit_seq = (sign(randn(num_bits,1)) + 1)/2;
s_k = bits_to_4QAM(bit_seq);

%% 3)
SNR_dB = 20; 
SNR_lin = 10^(SNR_dB/10);
Es = 2;               
sigma2_n = Es /SNR_lin;   

% CN(0, sigma2_n)
n1 = sqrt(sigma2_n/2) * (randn(N,1) + 1i*randn(N,1));
n2 = sqrt(sigma2_n/2) * (randn(N,1) + 1i*randn(N,1));

% outputs
r1 = h1_block .* s_k + n1;
r2 = h2_block .* s_k + n2;

%% 4)
% R = (h^H/h_norm)Y 
hHY = conj(h1) * r1 + conj(h2) * r2; 
h_norm = sqrt(abs(h1)^2 + abs(h2)^2);

R = hHY / h_norm;
s_hat = sign(real(R)) + 1i * sign(imag(R));
ser = sum(s_hat ~= s_k) / length(s_k);

%% 5)
SNR_dB = 0:2:20;
K = 2000; 
N = 200;

BER_MRC_sim = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    
    total_bit_errors = 0;
    total_bits = 0;
    
    for k = 1:K
        % 1) 
        h1 = (randn(1,1) + 1i*randn(1,1)) / sqrt(2);
        h2 = (randn(1,1) + 1i*randn(1,1)) / sqrt(2);
      
        h1_block = h1*ones(N,1);
        h2_block = h2*ones(N,1);
        
        % 2) 
        num_bits = 2*N;
        bit_seq = (sign(randn(num_bits,1)) + 1)/2;
        s_k = bits_to_4QAM(bit_seq);
        
        % 3) 
        Es = 2;                            
        SNR_lin = 10^(SNR_dB(i)/10);           
        sigma2_n = Es / SNR_lin;               
        
        n1 = sqrt(sigma2_n/2) * (randn(N,1) + 1i*randn(N,1));
        n2 = sqrt(sigma2_n/2) * (randn(N,1) + 1i*randn(N,1));
        
        r1 = h1_block .* s_k + n1;
        r2 = h2_block .* s_k + n2;
        
        % 4)
        hHY = conj(h1) * r1 + conj(h2) * r2;
        h_norm = sqrt(abs(h1)^2 + abs(h2)^2);
        R = hHY / h_norm;
        
        % 5) 
        s_hat = sign(real(R)) + 1i * sign(imag(R));
        
        bit_seq_hat = zeros(num_bits,1);
        bit_seq_hat(1:2:end) = (real(s_hat) < 0);
        bit_seq_hat(2:2:end) = (imag(s_hat) < 0);
    
        bit_errors = sum(bit_seq ~= bit_seq_hat);
        total_bit_errors = total_bit_errors + bit_errors;
        total_bits = total_bits + num_bits;
        
    end
    
    BER_MRC_sim(i) = total_bit_errors / total_bits;
    
end

%% 6) 
SNR_lin = 10.^(SNR_dB/10);      

m = sqrt( SNR_lin ./ (2 + SNR_lin) );
BER_MRC_theory = ((1 - m)/2).^2 .* (2 + m);

BER_MRC_highSNR = 3 ./ (4 * (SNR_lin.^2));

% plot
figure;
semilogy(SNR_dB, BER_MRC_sim,     'bo-','LineWidth',1.5); hold on;
semilogy(SNR_dB, BER_MRC_theory,  'k--','LineWidth',2);
semilogy(SNR_dB, BER_MRC_highSNR, 'r','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR (4-QAM, MRC)');
legend('Simulated MRC','Theoretical MRC','High-SNR approx', ...
       'Location','southwest');

%% 7)
BER_TB_sim = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    
    total_bit_errors = 0;
    total_bits = 0;
    
    for k = 1:K
        % 1) 
        h1 = (randn(1,1) + 1i*randn(1,1)) / sqrt(2);
        h2 = (randn(1,1) + 1i*randn(1,1)) / sqrt(2);
        h_norm = sqrt(abs(h1)^2 + abs(h2)^2);
        
        % 2) 
        num_bits = 2*N;
        bit_seq = (sign(randn(num_bits,1)) + 1)/2;
        s_k = bits_to_4QAM(bit_seq);
        
        % 3) 
        Es = 2;                          
        SNR_lin = 10^(SNR_dB(i)/10);          
        sigma2_n = Es / SNR_lin;
        
        n = sqrt(sigma2_n/2) * (randn(N,1) + 1i*randn(N,1));
        
        % 4) transmit beamforming
        r = h_norm * s_k + n;
        
        % 5)
        s_hat = sign(real(r)) + 1i * sign(imag(r));
        
        bit_seq_hat = zeros(num_bits,1);
        bit_seq_hat(1:2:end) = (real(s_hat) < 0);
        bit_seq_hat(2:2:end) = (imag(s_hat) < 0);
        
        bit_errors = sum(bit_seq ~= bit_seq_hat);
        total_bit_errors = total_bit_errors + bit_errors;
        total_bits = total_bits + num_bits;
    end
    
    BER_TB_sim(i) = total_bit_errors / total_bits;
end

%% 8)
figure;
semilogy(SNR_dB, BER_MRC_sim, 'bo-','LineWidth',1.5,'MarkerSize',5); hold on;
semilogy(SNR_dB, BER_TB_sim,  'mo-','LineWidth',1.5,'MarkerSize',5);
semilogy(SNR_dB, BER_MRC_theory, 'k--','LineWidth',2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR: 1x2 MRC vs 2x1 Transmit Beamforming (4-QAM)');
legend('1x2 MRC (sim)','2x1 TB (sim)','Theory (order 2)', ...
       'Location','southwest');
   
%% 9)
BER_Alamouti_sim = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    
    total_bit_errors = 0;
    total_bits = 0;
    
    for k = 1:K
        % 1) 
        h1 = (randn(1,1) + 1i*randn(1,1)) / sqrt(2);
        h2 = (randn(1,1) + 1i*randn(1,1)) / sqrt(2);
        
        % 2)
        num_bits = 2*N;
        bit_seq = (sign(randn(num_bits,1)) + 1)/2;
        s_k = bits_to_4QAM(bit_seq); % E{|s_k|^2} = Es
        
        s1 = s_k(1:2:end) / sqrt(2); % E{|s1|^2} = Es/2
        s2 = s_k(2:2:end) / sqrt(2); % E{|s2|^2} = Es/2
        L  = length(s1); % num of Alamouti blocks
        
        % 3) 
        Es = 2;                       
        SNR_lin = 10^(SNR_dB(i)/10);       
        sigma2_n = Es / SNR_lin;           
        
        n1 = sqrt(sigma2_n/2) * (randn(L,1) + 1i*randn(L,1));  
        n2 = sqrt(sigma2_n/2) * (randn(L,1) + 1i*randn(L,1));  
        
        % 4) Alamouti code
        s1_star = conj(s1);
        s2_star = conj(s2);
        
        y1 = h1 * s1 + h2 * s2 + n1;
        y2 = h1 * (-s2_star) + h2 * s1_star + n2;
        
        h_norm = sqrt(abs(h1)^2 + abs(h2)^2);
        R1 = (conj(h1).*y1 + h2.*conj(y2)) / h_norm;   
        R2 = (conj(h2).*y1 - h1.*conj(y2)) / h_norm; 
        
        % 5) 
        s1_hat = sign(real(R1)) + 1i * sign(imag(R1));
        s2_hat = sign(real(R2)) + 1i * sign(imag(R2));
        
        s_hat = zeros(N,1);
        s_hat(1:2:end) = s1_hat;
        s_hat(2:2:end) = s2_hat;
        
        bit_seq_hat = zeros(num_bits,1);
        bit_seq_hat(1:2:end) = (real(s_hat) < 0);
        bit_seq_hat(2:2:end) = (imag(s_hat) < 0);
        
        bit_errors       = sum(bit_seq ~= bit_seq_hat);
        total_bit_errors = total_bit_errors + bit_errors;
        total_bits       = total_bits + num_bits;
        
    end
    
    BER_Alamouti_sim(i) = total_bit_errors / total_bits;
end

%% 10) 
figure;
semilogy(SNR_dB, BER_MRC_sim,      'bo-','LineWidth',1.5,'MarkerSize',5); hold on;
semilogy(SNR_dB, BER_TB_sim,       'mo-','LineWidth',1.5,'MarkerSize',5);
semilogy(SNR_dB, BER_Alamouti_sim, 'go-','LineWidth',1.5,'MarkerSize',5);
semilogy(SNR_dB, BER_MRC_theory,   'k--','LineWidth',2);

grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR: 1x2 MRC, 2x1 TB, 2x1 Alamouti (4-QAM, Rayleigh block fading)');
legend('1x2 MRC (sim)', '2x1 TB (sim)', '2x1 Alamouti (sim)', ...
       'Theory (order 2)', 'Location','southwest');
   
%=====================================================================================
%% 2ND PART 

N = 200; 
SNR_dB = 0:2:20;
K = 2000;              

constellation = [ 1+1i; 1-1i; -1+1i; -1-1i ]; % 4-QAM

BER_ML = zeros(size(SNR_dB));
BER_DECOR = zeros(size(SNR_dB));

for i_snr = 1:length(SNR_dB)
                    
    total_bit_errors_ML = 0;
    total_bit_errors_DECOR = 0;
    total_bits = 0;
    
    for k = 1:K   
        % 1) 
        h11 = (randn + 1i*randn)/sqrt(2);   % Tx1 -> Rx1
        h12 = (randn + 1i*randn)/sqrt(2);   % Tx2 -> Rx1
        h21 = (randn + 1i*randn)/sqrt(2);   % Tx1 -> Rx2
        h22 = (randn + 1i*randn)/sqrt(2);   % Tx2 -> Rx2
        
        % 2) 
        num_bits = 2*N;                      
        bits1 = (sign(randn(num_bits,1)) + 1)/2;
        bits2 = (sign(randn(num_bits,1)) + 1)/2;
        
        X1 = bits_to_4QAM(bits1);           
        X2 = bits_to_4QAM(bits2);            
        
        % 3) 
        Es = 4;                             
        SNR_lin = 10^(SNR_dB(i_snr)/10);
        sigma2_n = Es / SNR_lin; 
    
        W1 = sqrt(sigma2_n/2) * (randn(N,1) + 1i*randn(N,1));
        W2 = sqrt(sigma2_n/2) * (randn(N,1) + 1i*randn(N,1));
        
        Y1 = h11*X1 + h12*X2 + W1;           
        Y2 = h21*X1 + h22*X2 + W2;          
        
        %% 4a) ML 
        H = [h11 h12;
             h21 h22];          
        
        Y = [Y1 Y2].';          
        
        [X1_hat_ML, X2_hat_ML] = ML_2x2_4QAM(H, Y, constellation);
        
        %% 4b) Decorrelator
        [X1_hat_DECOR, X2_hat_DECOR] = decorrelator_2x2_4QAM(H, Y);
        
        %% 5)
        % --- ML ---
        bits1_hat_ML = zeros(num_bits,1);
        bits2_hat_ML = zeros(num_bits,1);
        
        bits1_hat_ML(1:2:end) = (real(X1_hat_ML) < 0);
        bits1_hat_ML(2:2:end) = (imag(X1_hat_ML) < 0);
        
        bits2_hat_ML(1:2:end) = (real(X2_hat_ML) < 0);
        bits2_hat_ML(2:2:end) = (imag(X2_hat_ML) < 0);
        
        % --- Decorrelator ---
        bits1_hat_DECOR = zeros(num_bits,1);
        bits2_hat_DECOR = zeros(num_bits,1);
        
        bits1_hat_DECOR(1:2:end) = (real(X1_hat_DECOR) < 0);
        bits1_hat_DECOR(2:2:end) = (imag(X1_hat_DECOR) < 0);
        
        bits2_hat_DECOR(1:2:end) = (real(X2_hat_DECOR) < 0);
        bits2_hat_DECOR(2:2:end) = (imag(X2_hat_DECOR) < 0);
        
        % Count errors
        err_ML    = sum(bits1 ~= bits1_hat_ML)    + sum(bits2 ~= bits2_hat_ML);
        err_DECOR = sum(bits1 ~= bits1_hat_DECOR) + sum(bits2 ~= bits2_hat_DECOR);
        
        total_bit_errors_ML    = total_bit_errors_ML    + err_ML;
        total_bit_errors_DECOR = total_bit_errors_DECOR + err_DECOR;
        
        total_bits = total_bits + 2*num_bits; % 2 streams × num_bits
    end
    
    BER_ML(i_snr)    = total_bit_errors_ML    / total_bits;
    BER_DECOR(i_snr) = total_bit_errors_DECOR / total_bits;
end

SNR_lin = 10.^(SNR_dB/10);
slope1 = 1 ./ SNR_lin;         
slope2 = 5 ./ (SNR_lin.^2);     

figure;
semilogy(SNR_dB, BER_ML,    'bo-','LineWidth',1.5,'MarkerSize',5); hold on;
semilogy(SNR_dB, BER_DECOR, 'mo-','LineWidth',1.5,'MarkerSize',5);
semilogy(SNR_dB, slope1,    'g-','LineWidth',1.5); % 1/SNR 
semilogy(SNR_dB, slope2,    'c-','LineWidth',1.5); % 5/SNR^2
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR for 2x2 system (4-QAM): ML vs Decorrelator');
legend('ML detector','Decorrelator','1/SNR','5/SNR^2', ...
       'Location','southwest');
