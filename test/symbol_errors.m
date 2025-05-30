function num_of_symbol_errors = symbol_errors(est_X, X)

num_of_symbol_errors=0;
for i=1:length(X)
    if est_X(i)~=X(i)
        num_of_symbol_errors=num_of_symbol_errors+1;
    end
end

end