%% Alexander Faust
%
% ECE 300 - Homework 2 - Decision Theory
%
% Spetember 25, 2023
%
clc; clear; close all;
%% MAP vs ML
% Part a)
b1 = 10;
r1 = 20;
b2 = 15;
r2 = 25;
% Refer to function segment at the bottom to see how prior, likelihood
% function, and a-posterior distribution is calculated

% Part d) - Test cases
%% Case 1:
case1_b1 = 1;
case1_r1 = 9;
case1_b2 = 1;
case1_r2 = 9;

% Generate decision vector and associated probability of error
case1_D_MAP = MAP(case1_b1, case1_r1, case1_b2, case1_r2);
case1_D_ML = ML(case1_b1, case1_r1, case1_b2, case1_r2);

case1_Perr_MAP = Perr(case1_D_MAP, case1_b1, case1_r1, case1_b2, case1_r2);
case1_Perr_ML = Perr(case1_D_ML, case1_b1, case1_r1, case1_b2, case1_r2);

%% Case 2:
case2_b1 = 4;
case2_r1 = 6;
case2_b2 = 1;
case2_r2 = 9;

% Generate decision vector and associated probability of error
case2_D_MAP = MAP(case2_b1, case2_r1, case2_b2, case2_r2);
case2_D_ML = ML(case2_b1, case2_r1, case2_b2, case2_r2);

case2_Perr_MAP = Perr(case2_D_MAP, case2_b1, case2_r1, case2_b2, case2_r2);
case2_Perr_ML = Perr(case2_D_ML, case2_b1, case2_r1, case2_b2, case2_r2);

%% Case 3:
case3_b1 = 1;
case3_r1 = 9;
case3_b2 = 4;
case3_r2 = 6;

% Generate decision vector and associated probability of error
case3_D_MAP = MAP(case3_b1, case3_r1, case3_b2, case3_r2);
case3_D_ML = ML(case3_b1, case3_r1, case3_b2, case3_r2);

case3_Perr_MAP = Perr(case3_D_MAP, case3_b1, case3_r1, case3_b2, case3_r2);
case3_Perr_ML = Perr(case3_D_ML, case3_b1, case3_r1, case3_b2, case3_r2);

% Display results:
P_err_MAP = [case1_Perr_MAP ; case2_Perr_MAP ; case3_Perr_MAP];
P_err_ML = [case1_Perr_ML ; case2_Perr_ML ; case3_Perr_ML];

D_MAP = [case1_D_MAP;case2_D_MAP;case3_D_MAP];
D_ML = [case1_D_ML;case2_D_ML;case3_D_ML];

Cases = ["Case I" "Case II" "Case III"]';

partd_table = table(Cases, D_MAP, P_err_MAP, D_ML, P_err_ML);
disp(partd_table);

% part e)
disp("Based on the results from the table above, it is clear that the decision rule for ML" + newline + ...
     "did not depend on the distribution of balls in Urn I. Since the rule was [1, 2] across all" + newline + ...
     "cases, this implies that when B2 occurs we CHOOSE B1 and for R2 occuring we CHOOSE R1. This" + newline + ...
     "means that when we pick R2 or B2 in Urn II we choose the same corresponding ball in Urn I,");

% part f)
disp("Furthermore, for the decision rule for MAP it is clear that the rule CAN depend on the" + newline + ...
     "distribution of balls in Urn I since the decision rule is NOT the same in all three cases." + newline);

% part g)
% List the red and blue balls in urn I & II for case 1:
urn1_case1 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 1];
urn2_case1 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 1];

% List the red and blue balls in urn I & II for case 2:
urn1_case2 = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2];
urn2_case2 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 1];

% List the red and blue balls in urn I & II for case 3:
urn1_case3 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 1];
urn2_case3 = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2];

% Run the simulation for each of the decision vectors determined by each
% algorithm:
error_case1_MAP = simulate(case1_D_MAP, urn1_case1, urn2_case1);
error_case1_ML = simulate(case1_D_ML, urn1_case1, urn2_case1);
error_case2_MAP = simulate(case2_D_MAP, urn1_case2, urn2_case2);
error_case2_ML = simulate(case2_D_ML, urn1_case2, urn2_case2);
error_case3_MAP = simulate(case3_D_MAP, urn1_case3, urn2_case3);
error_case3_ML = simulate(case3_D_ML, urn1_case3, urn2_case3);

% Store error estimates in a vector to display in the table:
Error_estimate_MAP = [error_case1_MAP ; error_case2_MAP ; error_case3_MAP];
Error_estimate_ML = [error_case1_ML ; error_case2_ML ; error_case3_ML];

simulation_results = table(Cases, P_err_MAP, Error_estimate_MAP, P_err_ML, Error_estimate_ML);
% Display results of the 10^5 simulation using the table:
disp(simulation_results);


%% Functions Createed
% Function to compute the prior probability of an event:
function pi_m = prior(b1, r1)
    % Compute prior probability as probability of picking a ball:
    pi_m = 1/(b1 + r1) .* [b1 r1];
end

% Function to compute the likelihood function of the event:
function P_rs = likelihood(b1, r1, b2, r2)
    % Compute probability that given b1 or r2 is "sent", what is the
    % probability of event b2 or r2 happening:
    P_rs = 1/(b2 + r2 + 1) .* [(b2 + 1) b2 ; r2 (r2+1)];
end

% Function to compute the a-posterior probability (using likelihood and
% prior:
function P_sr = posterior(b1, r1, b2, r2)
    % Compute a-posteriori:
    P_rs = likelihood(b1, r1, b2, r2);
    pi_m1 = [prior(b1, r1) ; prior(b1, r1)];
    pi_m2 = [prior(b2, b2) ; prior(r2, r2)];
    P_sr = P_rs .* pi_m1 ./ pi_m2;
end

% Function to compute the decision vector for the MAP algorithm:
function D_map = MAP(b1, r1, b2, r2)
    % Obtain Decision if Blue from urn 2:
    vec1 = likelihood(b1, r1, b2, r2);
    
    % Obtain the index where max occurs (1 => "blue" 2 => "red")
    blue_check = vec1(1, :) .* prior(b1, r1);
    if blue_check(1) < blue_check(2)
        check1 = 2;
    else
        if blue_check(1) == blue_check(2)
            % If equal, 1 or 2, doesnt matter
            check1 = randsample(2, 1);
        else
            check1 = 1;
        end
    end

    red_check = vec1(2, :) .* prior(b1, r1);
    if red_check(1) < red_check(2)
        check2 = 2;
    else
        if red_check(1) == red_check(2)
            % If equal, 1 or 2, doesnt matter
            check2 = randsample(2, 1);
        else
            check2 = 1;
        end
    end
    
    % Assign decision vector accordingly:
    D_map = [check1 check2];
end

% Function to compute the decision vector for the ML algorithm:
function D_ml = ML(b1, r1, b2, r2)
    % Obtain Decision if Blue from urn 2:
    vec1 = likelihood(b1, r1, b2, r2);
    
    % Obtain the index where max occurs (1 => "blue" 2 => "red")
    blue_check = vec1(1, :);
    if blue_check(1) < blue_check(2)
        check1 = 2;
    else
        if blue_check(1) == blue_check(2)
            % If equal, 1 or 2, doesnt matter
            check1 = randsample(2, 1);
        else
            check1 = 1;
        end
    end

    red_check = vec1(2, :);
    if red_check(1) < red_check(2)
        check2 = 2;
    else
        if red_check(1) == red_check(2)
            % If equal, 1 or 2, doesnt matter
            check2 = randsample(2, 1);
        else
            check2 = 1;
        end
    end
    
    % Assign decision vector accordingly:
    D_ml = [check1 check2];
end

% Function to compute the probaility of error for each possible decision
% vector based on logic:
function P_err = Perr(v, b1, r1, b2, r2)
    P_rs = likelihood(b1, r1, b2, r2);
    pi_br = prior(b1, r1);
    if v(1) == 1 % => we have B2 | B1
        if v(2) == 1 % => we have R2 | B1
            % Add (~ | R1)'s to P_err 
            P_err = P_rs(1,2) * pi_br(2) + P_rs(2, 2) * pi_br(2);
        elseif v(2) == 2 % => we have R2 | R1
            % Add P(B2|R1) + P(R2|B1)
            P_err = P_rs(1, 2) * pi_br(2) + P_rs(2, 1) * pi_br(1);
        end
    elseif v(1) == 2 % => we have B2 | R1
        if v(2) == 1 % => we have R2 | B1
            % Add P(R2|R1) + P(B2|B1)
            P_err = P_rs(2, 2) * pi_br(2) + P_rs(1, 1) * pi_br(1);
        elseif v(2) == 2 % => we have R2 | R1
            % Add P(R2|B1) + P(B2|B1)
            P_err = P_rs(2, 1) * pi_br(1) + P_rs(1, 1) * pi_br(1);
        end
    end
end

% Function to compute an estimated error based on the decision vector:
function est_error = simulate(D, urn1, urn2)
% Run the experiment 10^5 times:
N = 1e5;

% Initialize error count to zero:
error_count = 0;
for i = 1:N
    % Create "fresh" urn to sample from on each iteration:
    urn1_ = urn1;
    urn2_ = urn2;
    % Choose a random index from urn1 to choose from:
    rand_index_urn1 = randi(size(urn1_));

    % Obtain the ball from the randomly chosen index:
    rand_ball_urn1 = urn1_(rand_index_urn1);
    
    % Add the ball to urn2:
    urn2_ = [urn2_, rand_ball_urn1];
    
    % Choose a random index from urn2 ro choose from:
    rand_index_urn2 = randi(size(urn2_));
    % Pick the ball at that randomly chosen location:
    rand_ball_urn2 = urn2_(rand_index_urn2);

    % Check the decision vector corresponding to this choice of ball in urn
    % 2. If "red" check position 2 in decision vector. If "blue" check
    % position 1 in decision vector. IF the RECORDED ball from urn 1 does
    % not match the decision, then increment the probability of error count
    if rand_ball_urn2 == 1
        % Check position 1 of decision vector
        if D(1) ~= rand_ball_urn1           
            error_count = error_count + 1;
        end
    else
        % Check position 2 of decision vector
        if D(2) ~= rand_ball_urn1
            error_count = error_count + 1;
        end
    end
    % Loop back to the beginning now that we have checked whether or not
    % the decision matches the actual ball picked
end
% Compute estimated error as the total errors observed divided by total
% number of simulations:
est_error = error_count/N;

end
