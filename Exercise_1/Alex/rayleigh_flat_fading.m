function h = rayleigh_flat_fading(N,b,Rh0)
% Inputs:
%   N  - Number of samples 
%   b  - AR coefficient (0 < b < 1)
%   Rh0 - Channel power E{|h[k]|^2} = 1
%
% Outputs:
%   h - Complex channel coefficients

    var_e = Rh0 * (1 - b^2);

    % initialize channel vector
    h = complex(zeros(N, 1));
    h(1) = 1;
    
    % e~CN(0,var_e)
    e = sqrt(var_e/2) * (randn(N, 1) + 1i * randn(N, 1));
    
    % generate AR-1 process
    for k = 2:N
        h(k) = b * h(k-1) + e(k);
    end
end
