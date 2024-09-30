%% Alexander Faust
% ECE 300 Problem Set 1

clc; close all; clear;

% Problem 1a:
% List k vector for M-ary QAM:
K = [2 4 6 8];
M = 2.^K(1:length(K));

% Generate M-ary QAM constellation plots for K vector:
generateQAM(K);

% Problem 1b:
% First calculate energy per symbol for each scheme:
QAM4_E_s = 1/M(1) * sum((norm(QAM4).^2));
QAM16_E_s = 1/M(2) * sum((norm(QAM16).^2));
QAM64_E_s = 1/M(3) * sum((norm(QAM64).^2));
QAM256_E_s = 1/M(4) * sum((norm(QAM256).^2));

% Calculate energy per bit for each QAM scheme:
QAM4_E_b = 1/K(1) * QAM4_E_s;
QAM16_E_b = 1/K(2) * QAM16_E_s;
QAM64_E_b = 1/K(3) * QAM64_E_s;
QAM256_E_b = 1/K(4) * QAM256_E_s;

% Problem 1c:
Eb_vec = [QAM4_E_b QAM16_E_b QAM64_E_b QAM256_E_b];

figure;
subplot(2,1,1);
plot(K,Eb_vec, '-x', 'MarkerSize', 10, 'MarkerEdgeColor', 'r');
title("E_b versus k");
xlabel("K values");
ylabel("Energy per bit");

% Problem 1d:
N = 2;              % 2 dimensions for each case of QAM
SNR_PerBit = 2./(2.^K);

subplot(2,1,2);
plot(K,SNR_PerBit, '-x', 'MarkerSize', 10, 'MarkerEdgeColor', 'r');
title("Spectral Efficiency vs K");
xlabel("K values");
ylabel("Spectral Efficiency");

%% Problem 2
% Problem 2d:
% List parameters:
W = 1;
T = 0.25;
f_max = 100;
f_min = 0;
freqs = linspace(f_min, f_max, 10000);

A = 1/T;
% Generate H(f) and X(f) vectors:
X_f = A*T .* sinc(freqs*T) .* exp(-1j*pi*freqs*T);
H_f = exp(-log(2)/2 .* (freqs/W).^2);

% Calculate Y(f):
Y_f = H_f .* X_f;
Y_db = 20 * log10(abs(Y_f));            % Convert to dB scale

figure;
plot(freqs, Y_db);
title("Magnitude spectrum of Y(f)");
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");
ylim([-60 0]);

% Problem 2e:
location = find(Y_db < -50);
B_0 = freqs(location(1));
disp("Frequency at which the spectrum is suppressed by 50dB: " + B_0 + newline);

% Problem 2f - part 4:
% Create time interval:
t = linspace(0, 10*T, 1000);
sigma = sqrt(log(2)) / (2*pi*W);
y_Q = A*(qfunc((t-T)/sigma) - qfunc(t/sigma));

figure;
plot(t, y_Q);
ylabel("|y(t)| [dB]");
xlabel("time (s)");
title("|y(t)| in terms of Q-function");

% Problem 2f - part 5:
thresh = 0.1*max(abs(y_Q));
LT_thresh_indexes = find(abs(y_Q) < thresh);
T0_index = LT_thresh_indexes(1);
T0 = t(T0_index);
disp("Time T0 where signal stays below 10% of its peak value: " + T0 + newline);

% Provlem 2f - part 6:
B0T0 = B_0 * T0;

disp("The 'time-bandwidth product' B0T0 is: " + B0T0 + " which describes the reduction" + newline + ...
     "in occupied bandwidth for a given spectral confinement, which is a consequence to" + newline + ...
     "spreading the signal in time");


%% Functions Created

function generateQAM(K)
    for k = K
        M = 2^k;
        n = sqrt(M(:));
        
        % Generate coordinates with distance 1 between them, depending on
        % the value of k given
        x = linspace(-n/2 + 0.5, n/2 - 0.5, n);
        
        % Center values about zero
        x = x - mean(x);

        % Generate imaginary and real components of the points
        [Re, Im] = meshgrid(x,x);
        symbols = Re(:) + 1i * Im(:);

        % Store each sequence of points for each 2^k for calculations:
        var = "QAM" + M(:);
        assignin("base", var, symbols);
       
        % Plot Constellation
        figure;
        scatter(real(symbols), imag(symbols), 'x');
        grid on;
        xlabel("In-Phase");
        ylabel("Quadrature");
        axis equal;
        title(M(:) + "-QAM Constellation");        
    end
end

