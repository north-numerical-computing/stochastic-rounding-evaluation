% Test for the unit circle on p. 51 of
% The Princeton Companion to Applied Mathematics,
% Princeton University Press, 2015.


% Clear chop options and reset PRNG seed
clear all
rng(500)

% Steps of ODE integration
N = 512;
% Set up data sampling period
sampling_period = 10;

% Get the coordinates of the exact solution
h_dp = 2*pi/300;
for i=1:301
    coordinates_exact(i, :) = [cos((i-1)*h_dp), -sin((i-1)*h_dp)];
end
coordinates = [];

global precision

options.round = 1; % RN
% Run two cases: with RN and with SR.
for k = 1:2
    switch k
        case 1, options.format = 'fp16'; options.subnormal = 1; sr = 0;
        case 2, options.format = 'fp16'; options.subnormal = 1; sr = 1;
        precision = 10;
    end
    
    fprintf('k = %1.0f, prec = %s, subnormal = %1.0f\n',...
        k,options.format,options.subnormal)
    chop([],options)

    % Initial values
    u_fp = chop(1);
    v_fp = chop(0);
    
    % Timestep size
    x_fp = chop(0);
    h_fp = chop(2*pi/N, options);

    j=0;
    % Solve the unit circle
    for i=1:N+1
        % Update the coordinates
        if (mod(i, sampling_period) == 0)
            j = j+1;
            coordinates(k, j, :) = [u_fp, v_fp];
        end
        [u_fp, v_fp] = Euler(0, h_fp, u_fp, v_fp, sr, options); 
        x_fp = chop(x_fp + h_fp);
        options.round = 1;
        chop([],options)
    end
end

h = plot(...
       coordinates(1, :, 1), coordinates(1, :, 2),'-', ...
       coordinates(2, :, 1), coordinates(2, :, 2),'--', ...
       coordinates_exact(:, 1), coordinates_exact(:, 2),'--k');
xlabel('u')
ylabel('v')
grid
legend('fp16 RN','fp16 SR', ...
       'Exact', 'Position', [0.4 0.4 0.1 0.2])
set(h,'LineWidth',1.2)
pbaspect([1 1 1]) % Set the ratio of the axes
xlim([-1.5 2])
ylim([-1.5 2])

function [u, v] = Euler(accurate, h, u_prev, v_prev, SR, options);
  if (accurate)
      u = u_prev + h*v_prev;
      v = v_prev - h*u_prev;
  elseif (SR)
      u = srSum(u_prev, srMultFMA(h, v_prev, ...
          rand(), options), rand(), options);
      v = srSum(v_prev, -srMultFMA(h, u_prev, ...
          rand(), options), rand(), options);
  else
      u = chop(u_prev + chop(h * v_prev, options), options);
      v = chop(v_prev - chop(h * u_prev, options), options);
  end
end

% Augmented multiplication algorithm based on FMA
function [sigma, error] = TwoMultFMA(a, b, options);
    sigma = chop(a*b, options);
    error = chop(a*b-sigma, options); % NOTE: This is an FMA.
end

% Multiplication with stochastic rounding
function f = srMultFMA(a, b, R, options);
    global precision
    % Switch to RZ
    options.round = 4;
    chop([], options);
    [sigma, error] = TwoMultFMA(a, b, options);
    Es = floor(log2(abs(sigma)));
    P = sign(error)*R*2^(Es-precision);
    f = chop(chop(error+P)+sigma);
    % Switch back to RTN
    options.round = 1;
    chop([], options);
end

% Augmented summation algorithm
function [sum, error] = TwoSum(a, b, options);
    sum = chop(a + b, options);
    a_trunc = chop(sum - b);
    b_trunc = chop(sum - a_trunc);
    a_error = chop(a - a_trunc);
    b_error = chop(b - b_trunc);
    error = chop(a_error + b_error);
end

% Addition with stochastic rounding
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
