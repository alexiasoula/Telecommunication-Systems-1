function bit_seq = generate_bits(N)

    arguments
        N   (1,1)  {mustBePositive, mustBeInteger}
    end

    bit_seq = randi([0 1], 1, 4*N);   % 0/1, array 1×(4N)
end
