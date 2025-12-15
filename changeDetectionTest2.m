close all

% Free parameters of the problem 
n = 128; % number of (equidistant) data points on [0,1] 
SNR = 15; % SNR used to compute variance of the iid Gaussian noise added to the measurements
% gamma = 1e-2

%% Set up the model 

 %(v(1) <= x & x< v(2)) .*( -6/pi ) +...
f1 = @(v) @(x) (x <v(1)) .* (x + pi) +...
    (v(1) <= x & x< v(2)) .*(-1/2*sin(6*x)) +...
    (x >= v(2)) .*(-1/2*sin(6*2)); %sin(-x+pi)

% Data points and signal values 
data_points = linspace(-1, 1, n)'*pi; % equidistant data points 
jumpLoc = [ -1/2  1/2]*pi;

% Creates one f1 function
f = f1(jumpLoc);

% creates noisy function y
x = f(data_points); % function values at grid points 
F = eye(n);
y = F*x;
noise_variances = norm(x)^2/n/(10^(SNR/10));
noise = sqrt(noise_variances)*randn(size(F,1),1);
y = y + noise; % noisy data 
snr = 10*log10(norm(x)^2/n/noise_variances);,

%% Plotting noisy signal 

% plot noisy y
figure;
plot(data_points, x, 'LineWidth',2); hold on;
plot(data_points, y,'LineWidth',1.5); hold on;
legend("${f}$", "${y}$", 'fontsize', 16, 'Interpreter', 'latex');
title("True vs Noisy Signal");
xlabel("x (domain)");
ylabel("Signal value"); 

%% Get residual solution 
cvx_begin quiet
   variable xhat(n)
   minimize( norm(F*xhat-y) )
cvx_end

x_ls = xhat;

% set up residual operator
p = 3;

T = LocalDetector(n,p);
S = getConMatrix('loc',n,p,1/4);

Res = (T-S);

% get solution for residual operator 
alpha = 2^(3/2)*(noise_variances)/sqrt((sum((Res * x_ls).^2)) / n);

cvx_begin quiet
    variable xhat(n)
    minimize( norm(F*xhat-y) + alpha*norm(Res*xhat,1) )
cvx_end

x_R = xhat;

% plot the residual solution
figure;
plot(data_points, x,'k--', 'LineWidth',2);hold on;
plot(data_points, x_R,'linewidth',1.5); hold on;
legend("${x_R}$", 'fontsize', 16, 'Interpreter', 'latex');
title("Residual Solution");

%% TV operator

% regularization order
q=1;
P = TV_operator( n, q); % regularization operator (matrix of finite dif) 
% P{2} = TV_operator( n, q+2); % regularization operator 

alpha = 3;%2^(3/2)*(noise_variances)/sqrt(sum((P{1}*x_ls).^2)/n);

cvx_begin quiet
  variable xhat(n)
  minimize( norm(F*xhat-y) + alpha*norm(P*xhat,1) )
cvx_end

x_TV = xhat;

% plot the residual solution
figure;
plot(data_points, x,'k--', 'LineWidth',1.5); hold on;
plot(data_points, x_TV,'LineWidth',2); hold on;
legend("${x_TV}$", 'fontsize', 16, 'Interpreter', 'latex');
title("TV Solution");

%% Edge map for TV
TV_edges = T*x_TV;
TV_edge_abs = abs(TV_edges);

figure;
plot(TV_edges, 'LineWidth',2); hold on;
plot(TV_edge_abs, 'LineWidth',1.5); hold on;
legend("${TV_edges}$","${TV_edges_abs}$", 'fontsize', 16, 'Interpreter', 'latex');
title("TV edge map and Abs TV edge map");

thresh = 2;
beta = 0.2;

jumps = [];   % dynamic list okay here

for j = 1:n
    if TV_edge_abs > thresh
        jumps(end+1) = j;   % append index
    end
end

% Now we know how many jumps:
intervals = zeros(length(jumps)-1,2);

for k = 1:length(jumps)-1
    bounds = [1, jumps, n+1];               % note n+1 so last becomes [jm, n]
    intervals = [bounds(1:end-1)', bounds(2:end)'-1];
end

x_sol = zeros(size(y));                
n = numel(x_TV);

for k = 1:size(intervals,1)

    L = intervals(k,1);              
    R = intervals(k,2);              
    n_s = R - L + 1;

    edge_metric = sum(TV_edge_abs(L:R), 1) / n_s;

    if edge_metric > beta
        x_sol(L:R) = x_R(L:R);
    else
        x_sol(L:R) = x_TV(L:R);
    end

end

figure;
plot(data_points, x_sol); hold on;
legend("${x_sol}$", 'fontsize', 16, 'Interpreter', 'latex');
title("Combination");

%% Notes:
% cross validation, looking at the alpha param in the residual function
% 





