function X = bits_to_4QAM(b)
  
    % reshape to [N/2 x 2] matrix
    bits = reshape(b, 2, [])';  
    
    % map bits to amplitudes (4-QAM mapping)
    map = [1 -1];  % 0 -> 1, 1 -> -1
    
    % convert binary to indices
    idx_I = bits(:,1) + 1;  % For I (first bit)
    idx_Q = bits(:,2) + 1;  % For Q (second bit)
    
    I_amp = map(idx_I);  % Get I amplitudes
    Q_amp = map(idx_Q);  % Get Q amplitudes
    
    % make complex symbols
    X = (I_amp + 1i*Q_amp);
    X = X(:); % column vector output
end