% Tests for various summation algorithms.
% Includes summation of random values and harmonic series.

% Clear environment and reset PRNG seed
clear all
rng(300, 'mrg32k3a') % Used by chop
% Create another stream aligned with the one above
s1 = RandStream('mrg32k3a', 'seed', 300);
% Create another stream for summation
s2 = RandStream('mrg32k3a', 'seed', 500);

% Set up number of iterations
N = 100000;
% Set up data sampling period
sampling_period = 100;

global precision

errors = [];
% Calculate iteration numbers at which errors are sampled
indices = []; 
for i = 1:N
    if (mod(i, sampling_period) == 0)
        indices = [indices, i];
    end
end

% Set up working floating-point arithmetic
options.format = 'fp16';
precision = 10;
options.subnormal = 1;
options.round = 1;
chop([], options);
% A vector of random values for summing
random_values = chop(rand(s2, N, 1)*0.01-0.002, options);

sum_double_recursive = 0;
sum_reduced_recursive = 0;
sum_reduced_compensated = 0;
sum_reduced_doubly_compensated = 0;
sum_reduced_cascaded = 0;
sum_reduced_stochastic1 = 0;
sum_reduced_stochastic2 = 0;
sum_reduced_stochastic3 = 0;
n = 1;

% Run the sum in the increasing order
for i = 1:N
    % Choose what to add here
    addend = random_values(i); % Random value sum
    %addend = chop(1/i); % Harmonic sum
    
    sum_double_recursive = sum_double_recursive + addend;
    sum_reduced_recursive = chop(sum_reduced_recursive + addend);
    sum_reduced_compensated = ...
        compensated(sum_reduced_compensated, addend, options);
    sum_reduced_doubly_compensated = ...
        doubly_compensated( ...
        sum_reduced_doubly_compensated, addend, options);
    [sum_reduced_cascaded, sum_cascaded_final] = ...
        cascaded_sum(sum_reduced_cascaded, addend, options);
    % Note, we are generating a random number here in order to use the
    % same random number for both algorithms below.
    R = rand(s1);
    sum_reduced_stochastic1 = ...
        srSum(sum_reduced_stochastic1, addend, R, options);
    options.round = 5;
    sum_reduced_stochastic3 = chop( ...
        sum_reduced_stochastic3 + addend, options);
    options.round = 1;
    chop([], options);
    
    % Sample the sums
    if (mod(i, sampling_period) == 0)
        errors(1, n) = sum_double_recursive;
        errors(2, n) = sum_reduced_recursive;
        errors(3, n) = sum_reduced_compensated;
        errors(4, n) = sum_reduced_doubly_compensated;
        errors(5, n) = sum_cascaded_final;
        errors(6, n) = sum_reduced_stochastic1;
        errors(7, n) = sum_reduced_stochastic3;
        n = n + 1;
    end
end

% Plotting
h = plot(indices, errors(1, :), '-', ...
             indices, errors(2, :), '-', ...
             indices, errors(3, :), '-', ...
             indices, errors(4, :), '-', ...
             indices, errors(5, :), '-', ...
             indices, errors(6, :), '-', ...
             indices, errors(7, :), '-');
xlabel('terms')
ylabel('value of the sum')
grid
legend('fp64 recursive','fp16 resursve', ...
       'fp16 compensated', 'fp16 doubly-compensated', ...
       'fp16 cascaded', 'fp16 stochastic', ...
       'fp16 stochastic chop')
set(h,'LineWidth',1.5)

% Compensated summation (Kahan's)
function f = compensated(s, x, options);
    persistent e
    if isempty(e)
        e = 0;
    end
    z = s;
    y = chop(x + e, options);
    [s, e] = fastTwoSum(s, y, options);
    f = s;
end

% Doubly Compensated summation (Priest's)
function f = doubly_compensated(s, x, options);
    persistent e
    if isempty(e)
        e = 0;
    end
    [y, u] = fastTwoSum(e, x, options);
    [t, v] = fastTwoSum(s, y, options);
    error = chop(u + v);
    [f, e] = fastTwoSum(t, error, options);
end

% Cascaded summation (Rump, Ogita and Oishi)
function [sum, sum_final] = cascaded_sum(s, x, options);
    persistent e
    if isempty(e)
        e = 0;
    end
    [sum, error] = TwoSum(s, x, options);
    e = chop(e + error, options);
    sum_final = chop(sum + e, options);
end

% FastTwoSum
function [sum, error] = fastTwoSum(a, b, options);
    sum = chop(a + b, options);
    b_trunc = chop(sum - a);
    error = chop(b - b_trunc);
end

% TwoSum
function [sum, error] = TwoSum(a, b, options);
    sum = chop(a + b, options);
    a_trunc = chop(sum - b);
    b_trunc = chop(sum - a_trunc);
    a_error = chop(a - a_trunc);
    b_error = chop(b - b_trunc);
    error = chop(a_error + b_error);
end

% Stochastic rounding summation
function f = srSum(a, b, R, options);
    global precision
    % Switch on RN
    options.round = 1;
    [sum, error] = TwoSum(a, b, options);
    % Switch to RZ
    options.round = 4;
    Es = floor(log2(abs(chop(a+b, options))));
    P = sign(error)*R*2^(Es-precision);
    % Switch to RD or RU.
    if (error >= 0)
        options.round = 3;
    else
        options.round = 2;
    end
    chop([], options);
    f = chop(chop(error+P)+sum);
    % Switch back to RN
    options.round = 1;
    chop([], options);
end

% Stochastic rounding summation (faster version)
function f = fastSrSum(a, b, R, options);
    global precision
    % Switch to RZ
    options.round = 4;
    chop([], options);
    [sum, error] = TwoSum(a, b, options);
    Es = floor(log2(abs(sum)));
    P = sign(error)*R*2^(Es-precision);
    f = chop(chop(error+P)+sum);
    % Switch back to RN
    options.round = 1;
    chop([], options);
end
