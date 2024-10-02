%% Alexander Faust
%
% ECE 300 - Problem Set 3 - Information Theory
%
% October 8, 2023
clc; clear; close all;
%% Question 4 - Graph Channel Capacity as a Function of p_E
p_E = linspace(0, 1, 100);
C = 1 - p_E;

figure;
plot(p_E, C);
title("Plot of Channel Capacity vs P_{E}");
ylabel("Channel Capacity [bits]");
xlabel("p");

%% Question 5 - Water Filling Algorithm
% List variances to run the algorithm with:
variances = [4 3 0.8 0.1];

% Question 5 - part a)
% Find minimum D s.t. R <= 2
% Trial and error to find lambdas s.t. R is above 2 and R below 2
lamb1 = 0.8;
lamb2 = 1.2;
[D1, R1, ~, ~] = calcDandR(lamb1, variances);
[D2, R2, ~, ~] = calcDandR(lamb2, variances);

% Display straddle points R1 and R2:
disp(R1);
disp(R2);

% Find a value for lambda s.t. R between 1.99 and 2:
% Sweep lambdas in the range between lamb1 and lmb2
lambdas = linspace(lamb1, lamb2, 1000);
best_lambda = [];
R = zeros(1, 1000);
D = zeros(1, 1000);
for i = 1:length(lambdas)
    lambda_iter = lamb1 + (lambdas(i) - lamb1);
    [D(i), R(i), ~, ~] = calcDandR(lambda_iter, variances);
    if R(i) <= 2
        if R(i) > 1.99
            best_lambda = [best_lambda, lambda_iter];
       
        end
    end
    
end

disp(best_lambda);
% Question 5 - part c)
% Create the graph of the R(D) curve between lambda 1 and lambda 2:
plotRD(R, D, lambdas);

% Question 5 - part d)
disp("Interpreting the R(D) graph in part c), as lambda increases from left to" + newline + ...
     "right, the distortion increases while the rate decreases. This makes sense" + newline + ...
     "intuitively for instance with image compression. Better looking images have" + newline + ...
     "less distortion, while increasing the rate they can be transfered due to this.");

% Question 5 - part e)
% Choose the middle element from the "optimal lambdas" computed earlier:
optimal_lambda = best_lambda(round(end/2));
[D_opt, R_opt, D_i, R_i] = calcDandR(optimal_lambda, variances);
disp("D calculated from optimal lambda: " + D_opt);
disp("R calculated from optimal lambda: " + R_opt);
disp("R_i for each vector component: ");
disp(R_i);
disp("D_i for each vector component: ");
disp(D_i);
disp("The third and fourth components of R_i are zeros!");


%% Functions Created
% Create function to change distortion D and rate R based on choice of
% lambda:
function [D, R, Di, Ri] = calcDandR(lambda, var)
    
    % Assign distortion vector D_i based on choice of lambda compared with
    % each variance value:
    Di = zeros(1, length(var));
    Ri = zeros(1, length(var));
    for i = 1 : length(var)
        if var(i) >= lambda
            Di(i) = lambda;
        else
            Di(i) = var(i);
        end
        Ri(i) = 1/2*log2(var(i) / Di(i));
    end
    % Compute D from the D_i vector:
    D = sum(Di);
    % Compute R & R_i:
    R = sum(Ri);

end

function plotRD(R, D, lambdas)
    figure;
    hold on;
    plot(D, R);
    title("Plot of R(D) Curve and With Lambda between \lambda_{1} = " + lambdas(1) + " & \lambda_{2} = " + lambdas(end));
    ylabel("R - Rate [bits]");
    xlabel("D - Distortion (mean-square)");
    plot([D(1), D(end)], [R(1), R(end)]);
    % Superimpose plot of increasing lambda:
    
    legend("R(D) curve", "Line Segment (D(\lambda_{1}), R(\lambda_{1})) -> (D(\lambda_{2}), R(\lambda_{2}))");
    hold off;
end
