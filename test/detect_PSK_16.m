function [est_X, est_bit_seq] = detect_PSK_16(Y)
% Detect 16-PSK symbols using Gray code
%
% INPUT:
%   Y : 2xN matrix of received I/Q components
%
% OUTPUT:
%   est_X         : 2xN matrix, estimated I/Q symbols
%   est_bit_seq   : 1x(4*N) estimated bit sequence (row vector)

    % Define Gray code
    gray_symbols = [
        0 0 0 0;
        0 0 0 1;
        0 0 1 1;
        0 0 1 0;
        0 1 1 0;
        0 1 1 1;
        0 1 0 1;
        0 1 0 0;
        1 1 0 0;
        1 1 0 1;
        1 1 1 1;
        1 1 1 0;
        1 0 1 0;
        1 0 1 1;
        1 0 0 1;
        1 0 0 0
    ];

    N = size(Y, 2);  % number of received symbols
    neighbours = zeros(1,16);
    % initialize vectors
    est_X = zeros(2, N);
    est_bit_seq = zeros(1, 4*N);

    % reference angles
    ref_angles = zeros(16,1);
    for i = 0:15
        ref_angles(i+1) = 2*pi*i/16;  
    end

    for n = 1:N
        % received I and Q
        I = Y(1,n);
        Q = Y(2,n);

        % calculate angle and normalize
        theta = atan2(Q, I);
        if theta < 0
            theta = theta + 2*pi;
        end
        nearest_neighbour = inf;
        idx = 1;
        % closest neighbour
        for j = 1:16
            neighbours(j) = sqrt( (I - cos((j-1)*2*pi/16))^2 + (Q - sin((j-1)*2*pi/16))^2 );  % euclidean distance 
            if nearest_neighbour > neighbours(j)
                nearest_neighbour = neighbours(j);
                idx = j;
            end
        end

        % estimated I and Q
        est_X(:,n) = [cos(ref_angles(idx)); sin(ref_angles(idx))];

        % append bits
        est_bit_seq(4*(n-1)+1 : 4*n) = gray_symbols(idx, :);
    end
end
