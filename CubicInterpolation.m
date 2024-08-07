function point = CubicInterpolation(z_row, mu)
    % This function returns an interpolated value of a point 
    % at a distance mu of second input sample. An interpolant
    % is approximated by third degree polynomial.
    %
    % Input: 
    % z_row -- a vector of input signal samples used for interpolation 
    % a row of 4 elements: z(n-1), z(n), z(n+1), z(n+2)  
    % mu -- fractional delay, |mu| < 1
    %
    % Output:
    % point -- a value of an interpolated sample.
    
    mu_col = [mu^3; mu^2; mu; 1]; % a column vector containing all powers of fractional delay mu (up to 3)
    B = [-1/6 1/2 -1/3 0; 1/2 -1 -1/2 1; -1/2 1/2 1 0; 1/6 0 -1/6 0]; % fixed coefficients for all input samples
    point = z_row * B * mu_col;
end