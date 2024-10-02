%% Alexander Faust
%
% ECE-300 Problem Set 4
%
% November 6, 2023
clc; clear; close all;

%% Question 5 - part a)
N = 10^7; 
SNR = 10:2:30;

phase = zeros(N, length(SNR));      
envelope = zeros(N, length(SNR));    

A = 1;
for i = 1: length(SNR)
    % List variance based on equation
    var = sqrt(10^(-SNR(i)/10));

    % Pick n_I and n_Q iid with N(0, 1/2*var)
    n_I = var * randn(N,1) / sqrt(2);
    n_Q = 1j * (var * randn(N,1) / sqrt(2));
    % Create total noise vector:
    n = n_Q + n_I;
    
    PDF = A + n; 
    % Phase of A + n is phase:
    phase(:,i) = angle(PDF);
    % Magnitude of A + n is envelope:
    envelope(:,i) = abs(PDF); 
end 


SNR_Samples = [10, 20, 30];
for i = 1:3 

    % Find indices where SNR vector is 10, 20, 30:
    index = find(SNR == SNR_Samples(i));

    % plot envelope histrogram 
    figure; 
    histogram(envelope(:,index), 'Normalization', 'pdf');
    title(['Histogram for SNR = ' num2str(SNR(index))]);
    
    %plot phase histrogram
    figure; 
    hold on;
    histogram(phase(:,index), 'Normalization', 'pdf');
    title(['Histogram for SNR = ' num2str(SNR(index))]);

    % Effie Bluestone helped me here:
    yline = ylim; % add in the y lim
    plot([10*(pi/180), 10*pi/180], yline, 'b-');
    plot([-10*pi/180, -10*pi/180], yline, 'b-');
    hold off; 
end 

%% Part 2

P_TH = zeros(1, length(SNR)); % pre allocate space 

for i = 1:length(SNR)

    % Evaluate if P when phase error exceeds 10 degrees:
    P_TH(i) = sum(abs(phase(:,i)) > 10 * pi/180) / N;
end

% Plot P_TH in log based scale
figure;
plot(SNR, log10(P_TH));
xlabel('SNR (dB)');
ylabel('log10(PTH)');
title('SNR vs log10(PTH)');

%% Part 3

% Plot SNR with 2dB increase:
figure;
hold on;
plot(SNR(2:end), diff(log10(P_TH)));
xlabel('SNR (dB)');
ylabel('Diff of log10(PTH)');
title('Diff of log10(PTH)');
hold off; 
