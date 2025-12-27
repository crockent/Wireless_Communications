function [BER_ML, BER_ZF] = calculate_ber_4qam_SM_ML_ZF(N, K, SNR_dB)

    const = [-1-1j -1+1j 1-1j 1+1j];
    bits_const = [0 0; 0 1; 1 0; 1 1];

    BER_ML = zeros(size(SNR_dB));
    BER_ZF = zeros(size(SNR_dB));

    for iSNR = 1:length(SNR_dB)
        EsN0 = 10^(SNR_dB(iSNR)/10);
        N0   = 1/EsN0;          % Es=1

        errML = 0;
        errZF = 0;
        totBits = 0;

        for k = 1:K
            % 1) 2x2 Rayleigh block-fading κανάλι
            H = (randn(2,2)+1j*randn(2,2))/sqrt(2);

            % 2) Δύο 4-QAM ακολουθίες
            bits1 = randi([0 1], 2*N, 1);
            bits2 = randi([0 1], 2*N, 1);

            s1 = bits_to_4qam(bits1);
            s2 = bits_to_4qam(bits2);

            X = [s1.'; s2.'];        % 2 x N

            % 3) Έξοδοι
            noise = sqrt(N0/2)*(randn(2,N)+1j*randn(2,N));
            Y = H*X + noise;

            % 4α) ML detection (εξάντληση 4^2=16 συνδυασμών)
            [c1,c2] = ndgrid(const,const);
            cand = [c1(:).'; c2(:).'];   % 2 x 16
            HY = H*cand;                 % 2 x 16

            idxML = zeros(1,N);
            for n = 1:N
                diff   = Y(:,n) - HY;
                metric = sum(abs(diff).^2,1);
                [~,idxML(n)] = min(metric);
            end
            Xhat_ML = cand(:,idxML);     % 2 x N

            % 4β) Decorrelator (ZF)
            Xhat_ZF = H\Y;               % pinv(H)*Y

            % Demapping σε bits
            % ML
            bhat1_ML = qam4_to_bits(Xhat_ML(1,:).', const, bits_const);
            bhat2_ML = qam4_to_bits(Xhat_ML(2,:).', const, bits_const);
            errML = errML + sum(bhat1_ML ~= bits1) + sum(bhat2_ML ~= bits2);

            % ZF: hard decision στο κοντινότερο σημείο
            Xhard_ZF = zeros(2,N);
            for n = 1:N
                [~,i1] = min(abs(Xhat_ZF(1,n)-const).^2);
                [~,i2] = min(abs(Xhat_ZF(2,n)-const).^2);
                Xhard_ZF(1,n) = const(i1);
                Xhard_ZF(2,n) = const(i2);
            end
            bhat1_ZF = qam4_to_bits(Xhard_ZF(1,:).', const, bits_const);
            bhat2_ZF = qam4_to_bits(Xhard_ZF(2,:).', const, bits_const);
            errZF = errZF + sum(bhat1_ZF ~= bits1) + sum(bhat2_ZF ~= bits2);

            totBits = totBits + 4*N;     % 2 streams * 2 bits/symbol * N
        end

        BER_ML(iSNR) = errML/totBits;
        BER_ZF(iSNR) = errZF/totBits;
    end
end