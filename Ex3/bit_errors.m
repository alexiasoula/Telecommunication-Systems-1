function num_bit_errors = bit_errors(est_bit_seq, b)
    % Ensure inputs are vectors and have the same orientation
    est_bit_seq = est_bit_seq(:);
    b = b(:);

    % Calculate number of bit errors using biterr
    num_bit_errors = sum(est_bit_seq ~= b);
end