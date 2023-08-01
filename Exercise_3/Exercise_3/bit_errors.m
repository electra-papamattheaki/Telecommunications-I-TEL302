function num_of_bit_errors = bit_errors(est_bit_seq, b)

    num_of_bit_errors = sum(est_bit_seq ~= b);
   
end

