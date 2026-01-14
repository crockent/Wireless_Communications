function symbols = bits_to_4qam(bits)
% Check input length
if mod(length(bits), 2) ~= 0
error('Length of bits must be even');
end
% Group bits into pairs
bit_pairs = reshape(bits, 2, []).';
real_part = 2*bit_pairs(:,1) - 1; % maps 0->-1, 1->+1
imag_part = 2*bit_pairs(:,2) - 1; % maps 0->-1, 1->+1
% Form complex symbols
symbols = real_part + 1j*imag_part;
end