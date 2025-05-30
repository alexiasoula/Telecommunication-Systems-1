function num_of_bit_errors = bit_error(est_bit_seq, b)
num_of_bit_errors=0;
    for i=1:length(b)
        if  est_bit_seq(i) ~= b(i)
         num_of_bit_errors=num_of_bit_errors+1;
       end
    end                  

end