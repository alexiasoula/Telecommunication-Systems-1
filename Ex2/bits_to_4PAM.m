%Function for A.4
function [X] = bits_to_4PAM(n)
    b1 = (sign(randn(n, 1)) + 1)/2;
    b2 = (sign(randn(n, 1)) + 1)/2;

    b = b2 * 10 + b1;

    for i = 1 : n
        if b(i) == 00
            X(i) = 3;
        elseif b(i) == 01
            X(i) = 1;
        elseif b(i) == 11
            X(i) = -1;
        else
            X(i) = -3;
        end
    end
end