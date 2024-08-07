function z = MyLinearConvolution(x, y)
    % This function returns the linear convolution of two input vectors.
    % Direct and FFT-based calculation methods are implemented. Second
    % method is used under the assumption that the discrete Fourier
    % transform of convolution of two sequences is a product of DFTs of
    % this sequences.
    %
    % Input: 
    % Vectors x and y. 
    %
    % Output:
    % z -- vector elements of which are the convoltuion of input vectors.

    Lx = length(x);
    Ly = length(y);

    N = Lx + Ly - 1;
    bit_depth = length(dec2bin(N));
    M = 2^bit_depth;

    % Determine the calculation method
    if 9*N^2 - 2*N + 1 < M * (27 * bit_depth - 8)
        % Direct convolution computation
        z = zeros(1, N);
        y_pad = [y(end:-1:1) zeros(1, Lx - 1)];
        for k = 1:N
            y_shift = circshift(y_pad, k-1);
            z(k) = x * y_shift(Ly:end).';
        end
    else
        % Fast Fourier transform based
        FFT_size = 2^bit_depth;
        X = fft(x, FFT_size);
        Y = fft(y, FFT_size);
        Z = X .* Y;
        z_ext = ifft(Z);
        z = z_ext(1:N);
    end
end

