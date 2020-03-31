% Testbench for ODE solvers with different rounding routines.
% Reworked version of https://epubs.siam.org/doi/10.1137/19M1251308

% Clear chop options and reset PRNG seed
clear options
rng(1)

% Set up some ODE conditions
a = 0; b = 0.015625;
y0 = 1.0;

global precision

% Exact solution to the exponential decay ODE
yexact = exp(-0.015625/20)*y0;

% Range of timesteps
nrange = round(10.^linspace(1, 3, 16));
m = length(nrange);

% Standard double precision. 
for j = 1:m
    n = nrange(j);
    x_dp = a;
    h_dp = (b-a)/n;
    y_dp = y0;

    for i=1:n
         y_dp = Euler(1, h_dp, x_dp, y_dp, 0, []);
         %y_dp = Midpoint(1, h_dp, x_dp, y_dp, 0, []);
         %y_dp = Heun(1, h_dp, x_dp, y_dp, 0, []);
        x_dp = x_dp + h_dp;
    end

    efp(j, 1) = abs(y_dp - yexact);
end

options.round = 1; % RN
% All chop formats.
for k = 1:6
    switch k
        case 1, options.format = 'b'; precision = 7; ...
            options.subnormal = 1; sr = 0;
        case 2, options.format = 'b'; precision = 7; ...
            options.subnormal = 1; sr = 1;
        case 3, options.format = 'h'; precision = 10; ...
            options.subnormal = 1; sr = 0;
        case 4, options.format = 'h'; precision = 10; ...
            options.subnormal = 1; sr = 1;
        case 5, options.format = 's'; precision = 23; ...
            options.subnormal = 1; sr = 0;
        case 6, options.format = 's'; precision = 23; ...
            options.subnormal = 1; sr = 1;
    end
    
    fprintf('k = %1.0f, prec = %s, subnormal = %1.0f\n', ...
        k,options.format,options.subnormal)
    chop([],options)
    
    % NOTE: We do not want to chop these, or any other constants with SR.
    a = chop(a); b = chop(b); y0 = chop(y0);
    
    for j = 1:m
        n = nrange(j);
        x_fp = chop(a);
        h_fp = chop((b-a)/n);
        y_fp = chop(y0);
        
        for i=1:n
            y_fp = Euler(0, h_fp, x_fp, y_fp, sr, options);
            %y_fp = Midpoint(0, h_fp, x_fp, y_fp, sr, options);
            %y_fp = Heun(0, h_fp, x_fp, y_fp, sr, options);
            %x_fp = chop(x_fp + h_fp, options);
        end
        efp(j,k+1) = abs(y_fp - yexact);
    end
end

h = loglog(...
       nrange,efp(:,1),'d-', ...
       nrange,efp(:,2),'x--', ...
       nrange,efp(:,3),'*--', ...
       nrange,efp(:,4),'o--', ...
       nrange,efp(:,5),'s-.', ...
       nrange,efp(:,6),'o-', ...
       nrange,efp(:,7),'-');
xlabel('$n$','Interpreter','latex')
ylabel('Error')
grid
legend('fp64','bfloat16 RN','bfloat16 SR', ...
       'fp16 RN','fp16 SR', 'fp32 RN',...
       'fp32 SR','Position',[0.69 0.6 0.1 0.2])
set(h,'LineWidth',1.2)

function f = decay_ODE(accurate, x, y, options);
    if (accurate)
        f = -y/20;
    else
        f = -chop(y/20, options);
    end
end

% ODE solvers
function f = Euler(accurate, h, x, y, SR, options);
  if (accurate)
      f = y + h*decay_ODE(1, x, y, options);
  elseif (SR)
      f =srSum(y, srMultFMA(h, decay_ODE(0, x, y, options), ...
          rand(), options), rand(), options);
  else
      f = chop(y + ...
        chop(h*chop(decay_ODE(0, x,y, options), options), options), ...
            options);
  end
end

function f = Midpoint(accurate, h, x, y, SR, options);
  if (accurate)
      f = y + h*decay_ODE(1, x, y + h*1/2*decay_ODE(1, x, y, options), ...
          options);
  elseif (SR)
      temp = srMultFMA(h, srMultFMA(1/2, decay_ODE(0, x, y, options), ...
          rand(), options), rand(), options);
      temp = srSum(y, temp, rand(), options);
      f = srSum(y, srMultFMA(h, decay_ODE(0, x, temp, options), ...
          rand(), options), rand(), options);
  else
      f = chop(y + ...
        chop(h*decay_ODE(0, x, ...
            chop(y + ...
                chop(chop(1/2*h, options)*decay_ODE(0, x, y, options), ...
                options)), options), options), options);
  end
end

function f = Heun(accurate, h, x, y, SR, options);
  if (accurate)
      y_temp = y + h*decay_ODE(1, x, y, options);
      f = y + h/2 * (decay_ODE(1, x, y, options) + ...
          decay_ODE(1, x, y_temp, options));
  elseif(SR)
      y_temp = srSum(y, srMultFMA(h, decay_ODE(0, x, y, options), ...
          rand(), options), rand(), options);
      f = srSum(y, srMultFMA(h/2, srSum(decay_ODE(0, x, y, options), ...
          decay_ODE(0, x, y_temp, options), rand(), options), rand(), ...
          options), rand(), options);
  else
      y_temp = chop(y + chop(h*decay_ODE(0, x, y, options), options));
      f = chop(y + chop(chop(h/2, options) * ...
          chop(decay_ODE(0, x, y, options) + ...
          decay_ODE(0, x, y_temp, options), options)));
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
