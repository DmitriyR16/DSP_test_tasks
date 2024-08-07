function [U1, U2] = CreateStateMatrices(x, K1, K2, Q)
% This function shapes two matrices of the state of an input
% signal with given history lengths and the order of nolinearity
%
% Input:
% 
% x - vector of input signal
% K1 and K2 are history lengths
% Q: a nonlinearity of an odd order up to 2*Q-1 is used
%
% Output:
%
% U1 and U2 are state matrices

    u1 = zeros(1, K1*Q); % an initial state vector of the first system
    u2 = zeros(1, K2*Q); % an initial state vector of the second system
    U1 = [];
    U2 = [];
    q = 0:Q-1;
    
    for n = 1:numel(x)
        u1 = circshift(u1, Q);
        u2 = circshift(u2, Q);
        % u_n = [x(n) x(n)*x(n)' x(n)^4*x(n)'];
        u_n = x(n) * (x(n)*x(n)').^q;
        u1(1:Q) = u_n;
        u2(1:Q) = u_n;
        U1 = [U1; u1];
        U2 = [U2; u2];
    end

end