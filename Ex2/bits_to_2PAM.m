%Function for A.2
function [X] = bits_to_2PAM(n)
    b = (sign(randn(n, 1)) + 1)/2;
    for i = 1 : n
        if b(i) == 0
            X(i) = 1;
        else 
            X(i) = -1;
        end
    end
end