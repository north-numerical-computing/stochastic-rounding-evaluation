% Test summation algorithms in native precision with stochastic rounding.

% Clear environment and reset PRNG seed
clear all
rng(1)

% Set up number of iterations
N = 1000000;
% Set up data sampling period
sampling_period = 1000;

errors = [];

% Calculate iteration numbers at which errors are sampled
indices = []; 
for i = 1:N
    if (mod(i, sampling_period) == 0)
        indices = [indices, i];
    end
end

% A vector of random fp64 values for summing
random_values = rand(N, 1)*2^-65;

sum_double_recursive = 1;
sum_double_compensated = 1;
sum_double_cascaded = 1;
sum_double_stochastic1 = 1;
sum_double_stochastic2 = 1;
n = 1;

% Run the sum in the increasing order
for i = 1:N
    addend = random_values(i);
    sum_double_recursive = sum_double_recursive + addend;
    sum_double_compensated = compensated(sum_double_compensated, addend);
    [sum_double_cascaded, sum_cascaded_final] = ...
        cascaded_sum(sum_double_cascaded, addend);
    R = rand();
    sum_double_stochastic1 = sradd(sum_double_stochastic1, addend, R);
    
    % Sample the absolute errors
    if (mod(i, sampling_period) == 0)
        errors(1, n) = sum_double_recursive-1;
        errors(2, n) = sum_double_compensated-1;
        errors(3, n) = sum_cascaded_final-1;
        errors(4, n) = sum_double_stochastic1-1;
        n = n + 1;
    end
end

h = plot(indices, errors(1, :), '-', ...
         indices, errors(2, :), '-', ...
         indices, errors(3, :), '--', ...
         indices, errors(4, :), '-');
xlabel('Iterations i')
ylabel('s - 1')
grid
legend('fp64 recursive', ...
       'fp64 compensated', 'fp64 cascaded', ...
       'fp64 stochastic v.1')
set(h,'LineWidth',1.5)
         
% Compensated summation
function f = compensated(s, x);
    persistent e
    if isempty(e)
        e = 0;
    end
    z = s;
    y = x + e;
    s = z + y;
    e = (z-s)+y;
    f = s;
end

% Cascaded summation (Rump, Ogita and Oishi)
function [sum, sum_final] = cascaded_sum(s, x);
    persistent e
    if isempty(e)
        e = 0;
    end
    [sum, error] = twosum(s, x);
    e = e + error;
    sum_final = sum + e;
end
