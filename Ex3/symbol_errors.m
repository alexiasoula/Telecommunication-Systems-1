function num_errors = symbol_errors(est_X, X)
    % Ensure inputs are vectors and have the same orientation
    est_X = est_X(:);
    X = X(:);
    
    % Calculate symbol errors using symerr
    [num_errors, ~] = symerr(est_X, X);
end
