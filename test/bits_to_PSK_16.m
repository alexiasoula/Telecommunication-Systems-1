function X = bits_to_PSK_16(bit_seq)
% a function for 16-PSK mapping using Gray code
%
% INPUT  bit_seq :
%   vector, length must be multiple of 4
%
% OUTPUT X :
%   2xN matrix, row1=I, row2=Q

len = length(bit_seq);
num_of_symbols = len/4;

X = zeros(2, num_of_symbols);
XI = X(1,:);
XQ = X(2,:);

x_idx = 1;

M = 16;
theta = 2*pi*(0:M-1)/M;
I_ref = cos(theta);
Q_ref = sin(theta);

for i = 1 : 4 : len-1
    % 0000
    if ( bit_seq(i) == 0 && bit_seq(i+1) == 0 && bit_seq(i+2) == 0 && bit_seq (i+3) == 0 )
        XI(x_idx) = cos(2*pi*0 / M);
        XQ(x_idx) = sin(2*pi*0 / M);

    % 0001
    elseif ( bit_seq(i) == 0 && bit_seq(i+1) == 0 && bit_seq(i+2) == 0 && bit_seq (i+3) == 1 )
        XI(x_idx) = cos(2*pi*1/ M);
        XQ(x_idx) = sin(2*pi*1/ M);

    % 0011
    elseif ( bit_seq(i) == 0 && bit_seq(i+1) == 0 && bit_seq(i+2) == 1 && bit_seq (i+3) == 1 )
        XI(x_idx) = cos(2*pi*2 / M);
        XQ(x_idx) = sin(2*pi*2 / M);
    
    % 0010
    elseif ( bit_seq(i) == 0 && bit_seq(i+1) == 0 && bit_seq(i+2) == 1 && bit_seq (i+3) == 0 )
        XI(x_idx) = cos(2*pi*3 / M);
        XQ(x_idx) = sin(2*pi*3 / M);
     
    % 0110
    elseif ( bit_seq(i) == 0 && bit_seq(i+1) == 1 && bit_seq(i+2) == 1 && bit_seq (i+3) == 0 )
        XI(x_idx) = cos(2*pi*4 / M);
        XQ(x_idx) = sin(2*pi*4 / M);
      
    % 0111
    elseif ( bit_seq(i) == 0 && bit_seq(i+1) == 1 && bit_seq(i+2) == 1 && bit_seq (i+3) == 1 )
        XI(x_idx) = cos(2*pi*5 / M);
        XQ(x_idx) = sin(2*pi*5 / M);    

   
    % 0101
    elseif ( bit_seq(i) == 0 && bit_seq(i+1) == 1 && bit_seq(i+2) == 0 && bit_seq (i+3) == 1 )
        XI(x_idx) = cos(2*pi*6 / M);
        XQ(x_idx) = sin(2*pi*6 / M);
 
    % 0100
    elseif ( bit_seq(i) == 0 && bit_seq(i+1) == 1 && bit_seq(i+2) == 0 && bit_seq (i+3) == 0 )
        XI(x_idx) = cos(2*pi*7 / M);
        XQ(x_idx) = sin(2*pi*7 / M);


    % 1100
    elseif ( bit_seq(i) == 1 && bit_seq(i+1) == 1 && bit_seq(i+2) == 0 && bit_seq (i+3) == 0 )
        XI(x_idx) = cos(2*pi*8 / M);
        XQ(x_idx) = sin(2*pi*8 / M);


    % 1101
    elseif ( bit_seq(i) == 1 && bit_seq(i+1) == 1 && bit_seq(i+2) == 0 && bit_seq (i+3) == 1 )
        XI(x_idx) = cos(2*pi*9 / M);
        XQ(x_idx) = sin(2*pi*9 / M);

         
    % 1111
    elseif ( bit_seq(i) == 1 && bit_seq(i+1) == 1 && bit_seq(i+2) == 1 && bit_seq (i+3) == 1 )
        XI(x_idx) = cos(2*pi*10 / M);
        XQ(x_idx) = sin(2*pi*10 / M);


    % 1110
    elseif ( bit_seq(i) == 1 && bit_seq(i+1) == 1 && bit_seq(i+2) == 1 && bit_seq (i+3) == 0 )
        XI(x_idx) = cos(2*pi*11 / M);
        XQ(x_idx) = sin(2*pi*11 / M);

    % 1010
    elseif ( bit_seq(i) == 1 && bit_seq(i+1) == 0 && bit_seq(i+2) == 1 && bit_seq (i+3) == 0 )
        XI(x_idx) = cos(2*pi*12 / M);
        XQ(x_idx) = sin(2*pi*12 / M);

    % 1011
    elseif ( bit_seq(i) == 1 && bit_seq(i+1) == 0 && bit_seq(i+2) == 1 && bit_seq (i+3) == 1 )
        XI(x_idx) = cos(2*pi*13 / M);
        XQ(x_idx) = sin(2*pi*13 / M);

    % 1001
    elseif ( bit_seq(i) == 1 && bit_seq(i+1) == 0 && bit_seq(i+2) == 0 && bit_seq (i+3) == 1 )
        XI(x_idx) = cos(2*pi*14 / M);
        XQ(x_idx) = sin(2*pi*14 / M);
        
    % 1000
    elseif ( bit_seq(i) == 1 && bit_seq(i+1) == 0 && bit_seq(i+2) == 0 && bit_seq (i+3) == 0 )
        XI(x_idx) = cos(2*pi*15 / M);
        XQ(x_idx) = sin(2*pi*15 / M);
    end
    x_idx = x_idx + 1;
end

X = [XI; XQ];

end
