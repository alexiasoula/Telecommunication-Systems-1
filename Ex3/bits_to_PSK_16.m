function X = bits_to_PSK_16(bit_seq)

    if ~isvector(bit_seq) || ~all(ismember(bit_seq,[0 1]))
        error('bit_seq must be a vector containing only 0s and 1s.');
    end
    L = numel(bit_seq);
    if mod(L,4) ~= 0
        error('Length of bit_seq (%d) is not a multiple of 4.', L);
    end
    N = L/4;
    bit_seq = bit_seq(:).';

    % 4 bits -> b (0–15)
    bits_mat = reshape(bit_seq, 4, N);        % 4×N
    b = bits_mat(1,:)*8 + bits_mat(2,:)*4 + bits_mat(3,:)*2 + bits_mat(4,:);

    % Gray
    g = bitxor(b, bitshift(b, -1));

    % 16-PSK
    theta = 2*pi * g / 16;                    % 1×N
    X = [cos(theta);                          % 2×N
         sin(theta)];
end