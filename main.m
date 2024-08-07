clc;
clear all;
close all;

%% Task 1. Linear convolution of two vectors with the same length

N_all = [2, 3, 16];

for i = 1:length(N_all)
    N_i = N_all(i);
    x = rand(1, N_i);
    y = rand(1, N_i);
    my_conv = MyLinearConvolution(x, y);
    std_conv = conv(x, y);
    fprintf(['For N = %d maximal difference between standard convolution and ' ...
        'mine equals %.4f\n'], N_i, max(abs(std_conv - my_conv)));
end

%% Task 2. Signal resampling

f_s = 1e4;          % sampling frequency in Hz
L = 4;
M = 5;
f_res = L/M * f_s;  % new sampling frequency in Hz
t = 0:(1/f_s):1;    % sequence of time samples
uniform_signal = randn(size(t));    % random signal with normal distributuion
% signal = lowpass(uniform_signal, f_s/6, f_s);
[b,a] = butter(4, 1/6); % generating filter coefficients for a 2nd order lowpass filter
signal = filter(b, a, uniform_signal); % applying filter to random signal
signal_spec = fft(signal) / sqrt(length(signal));
nu_s = linspace(-1/2, 1/2, length(signal)); % normalised frequency range

% Plot origin signal's spectrum 
figure;
plot(nu_s, abs(fftshift(signal_spec)));
grid on;
title('Origin signal spectrum'); 
xlabel('\nu'); 
ylabel('|Z|');

iter_num = floor(length(signal) / M);
resampled_signal = [];
% In blocks of 5 samples of input signal replace 4 samples with 3 
% interpolated ones iteratively
for i = 0 : iter_num-1
    interp_point_1 = CubicInterpolation(signal(i*M+1:i*M+4), 1/4);
    interp_point_2 = CubicInterpolation(signal(i*M+2:i*M+5), 1/2);
    interp_point_3 = CubicInterpolation(signal(i*M+3:i*M+6), 3/4);
    resampled_signal = [resampled_signal signal(i*M+1) interp_point_1 interp_point_2 interp_point_3];
end
% Measuring a quality of resampling
resampled_signal_spec = fft(resampled_signal) / sqrt(length(resampled_signal));
nu_res = linspace(-1/2, 1/2, length(resampled_signal));
diff = signal_spec(1:2000) - resampled_signal_spec(1:2000);
NMSE = 10*log10(diff * diff' / (signal_spec(1:2000) * signal_spec(1:2000)'));
fprintf(['Normalised mean square error between spectrums of origin and' ...
        ' resampled signals equals %.4f dB \n'], NMSE);

% Plot spectrum of origin and resampled signals 
figure;
plot(nu_s, abs(fftshift(signal_spec)));
hold on;
plot(nu_res, abs(fftshift(resampled_signal_spec)));
grid on;
title('Spectrum of signals');  
xlabel('\nu'); 
ylabel('|Signal spectrum|');
legend('|Z|, F_s', '|Y|, F_{res}');


%% Task3. Adaptive filtering
Bit_in_sym = 4;                 
Mod_type = '16QAM';
N = 1e5 * Bit_in_sym; % number of bits to be transmitted
bits = randi([0 1], 1, N); % bit sequence generation
symbols = Modulation(bits, Mod_type).'; % signal at the input of an estimated system
train_symbols = symbols(1:1e4); % data block for training
test_symbols = symbols(1e4+1:end); % data block for validation 
Q = 3; % P = 2q+1 is an order of nonlinearity, q = 0, ..., Q-1
K_syst = 10; % the length of a nonlinear model's memory
SNR_all = -15:3:30; 
NMSE_all = [];
K_filt_all = [];
%% Main loop
for i_snr = 1:numel(SNR_all)
    SNR = SNR_all(i_snr);
    system_coeffs = (rand(K_syst*Q, 1) + rand(K_syst*Q, 1) * 1j) / sqrt(2); % coefficients of the system
    K_filt = 0; % a length of the filter's memory
    est_coeffs = zeros(Q, 1); % an initial vector of filter' coefficients
    norm_err_pow = 1; % an initial value of normalised error power
    norm_err_pow_diff = 1; % an initial value of difference between normalised error powers 
    % on two adjacent steps

    % Adaptive filter training loop 

    while norm_err_pow_diff > 0
        K_filt = K_filt + 1;
        norm_err_pow_prev = norm_err_pow;

        % Forming the matrix of state vectors and the signal at the output of the system

        [U_syst_train, U_filt_train] = CreateStateMatrices(train_symbols, K_syst, K_filt, Q);
        y_train = U_syst_train * system_coeffs;
        y_noised = Add_WGNoise(y_train, SNR);
        
        % Estimate system coefficients by LS method and compute power of
        % error vector
    
        est_coeffs_prev = est_coeffs;
        est_coeffs = inv(U_filt_train' * U_filt_train + eye(size(U_filt_train, 2))*1e-2) * U_filt_train' * y_noised;
        y_adapt = U_filt_train * est_coeffs;
        err = y_noised - y_adapt;
        norm_err_pow = err' * err / (y_noised'*y_noised);
        norm_err_pow_diff = norm_err_pow_prev - norm_err_pow;
    end

    K_filt = K_filt - 1;
    K_filt_all = [K_filt_all K_filt];
    filter_coeffs = est_coeffs_prev;
    
    % Check results of adaptive filtering
    
    [U_syst_test, U_filt_test] = CreateStateMatrices(test_symbols, K_syst, K_filt, Q);
    y_test = U_syst_test * system_coeffs;
    y_adapt = U_filt_test * filter_coeffs;
    err = y_test - y_adapt;
    NMSE = 10*log10(err' * err / (y_test'*y_test));
    NMSE_all = [NMSE_all NMSE];
end

%% Plot NMSE
plot(SNR_all, NMSE_all, '-r*', 'LineWidth',2);
grid on;
title('The length of nonlinear model memory equals ', num2str(K_syst)); 
xlabel('SNR, dB'); 
ylabel('NMSE, dB');