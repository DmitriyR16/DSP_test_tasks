function noised_signal = Add_WGNoise (signal, SNR)
% This function returns noised signal with given signal to noise ratio
%
% Input:
%
% signal - a vector of input signal samples
% SNR - signal to (additive white Gaussian) noise ratio
%
% Output:
%
% noised_signal - a vector of noised signal samples

    % Compute standard deviation
    sigma = sqrt(sum(signal.*conj(signal))/(length(signal)*2*10^(SNR/10)));
    % Add white Gaussian noise
    noised_signal = signal + (1 + 1j) * normrnd(0, sigma, size(signal)) / sqrt(2);
end
