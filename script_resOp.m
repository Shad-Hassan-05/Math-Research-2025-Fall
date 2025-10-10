%% script_denoising_1d
%
% Description: 
%  Script to reconstruct a piecewise constant signal based on noisy blurred MMV data. 
% 
% Author: Yao Xiao 
% Date: Sept 01, 2024 
% Edited: Shad Hassan 
% Date: Oct 10, 2025
%

clear; close all; clc; % clean up
%warning('off','all') % in case any of the warnings become too anoying 


%% Free parameters 

% Free parameters of the problem 
n = 128; % number of (equidistant) data points on [0,1] 
SNR = 5; % SNR used to compute variance of the iid Gaussian noise added to the measurements
% gamma = 1e-2
J = 1; % number of (different) signals
rng('default'); rng(1 ,'twister'); % to make the results reproducable

beta = 10
% input('input beta: ');



%% Set up the model 
fconst = @(v) @(x) (x <v(1)) .* (1.5) +...
    (v(1) <= x & x< v(2)) .*( -6/pi ) +...
    (x >= v(2)) .*( 1.5);
 %(v(1) <= x & x< v(2)) .*(-1/2*sin(6*x)) +...
f1 = @(v) @(x) (x <v(1)) .* (x + pi) +...
    (v(1) <= x & x< v(2)) .*( -6/pi ) +...
    (x >= v(2)) .* sin(-x+pi);

% Data points and signal values 
data_points = linspace(-1, 1, n)'*pi; % equidistant data points 
jumpLoc = [ -1/2  1/2]*pi;

% Creates one f1 function
fun{1} = f1(jumpLoc);


% Regularization operator
% creates and nxn matrix with +1 and -1 on diagonals
% 
q=1; % regularization order
P{1} = TV_operator( n, q); % regularization operator (matrix of finite dif) 
P{2} = TV_operator( n, q+2); % regularization operator 


% F_conv = construct_F_deconvolution(n,gamma);
% r = 0.1
% F_us = construct_F_undersampling(n,gamma);

% Adds noise to the clean data
figure;
tlo = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for j=1:J 
    x{j} = fun{j}(data_points); % function values at grid points 
    
    


    F{j} = eye(n);

    y{j} = F{j}*x{j};
    noise_variances{j} = norm(x{j})^2/n/(10^(SNR/10));
    noise{j} = sqrt(noise_variances{j})*randn(size(F{j},1),1);
    y{j} = y{j} + noise{j}; % noisy data 
    snr = 10*log10(norm(x{j})^2/n/noise_variances{j})

    nexttile;
    plot(data_points, x{j}, 'LineWidth',2); hold on;
    plot(data_points, y{j},'LineWidth',1.5); hold on;
    xlim([data_points(1), data_points(end)]);
    if j ==1 legend("${f}$", "${y}$", 'fontsize', 16, 'Interpreter', 'latex');end
    title("True vs Noisy Signal");
    xlabel("x (domain)");
    ylabel("Signal value");


end  




%% edge detectors
% Grid and frequency setup
% builds the residual operator. S is gloabal fourier, T is the local diff.
% this is what we compair to the TV method. 
p = 0;


T = LocalDetector(n,p);
S = LocalDetector(n,3);
% S = getConMatrix('loc',n,p,1/4);
S_1 = getConMatrix('loc',n,3,1/4);

Res = (T-S);

%% Plotting the space of LocalDetector
figure;
tlo = tiledlayout(J,3,'TileSpacing','compact','Padding','compact');

for j = 1:J
    
    I = eye(n)

    sigma2 = norm(x{j})^2 / n / (10^(SNR/10));

    B = (1 + sigma2)*I

    sol = inv(B' * B)*B'*y{j}

    solT = T * y{j};
    
    solS = S *  y{j};

    nexttile;
    plot(data_points, sol, 'LineWidth', 2); hold on;
    xlim([data_points(1), data_points(end)]);
    if j == 1
        legend("${Tb}$", "${y}$", 'fontsize', 16, 'Interpreter', 'latex');
    end
    title("Space of Local Edge Detector");
    xlabel("x (domain)");
    ylabel("Signal value");

    nexttile;
    plot(data_points, solS, 'LineWidth', 2); hold on;
    xlim([data_points(1), data_points(end)]);
    if j == 1
        legend("${Sb}$", "${y}$", 'fontsize', 16, 'Interpreter', 'latex');
    end
    title("Space of Global Edge Detector");
    xlabel("x (domain)");
    ylabel("Signal value");

end


%% reconstruction section
% three data reconstruction methods, first is normal least squares.
% second is the TV operator, and the last is our residual operator. 
figure;
tlo = tiledlayout(J,3,'TileSpacing','compact','Padding','compact');

for j = 1:J
    cvx_begin quiet
        variable xhat(n)
        minimize( norm(F{j}*xhat-y{j}) )
    cvx_end

    x_ls = xhat;

    alpha = 2^(3/2)*(noise_variances{j})/sqrt(sum((P{1}*x_ls).^2)/n);

    cvx_begin quiet
        variable xhat(n)
        minimize( norm(F{j}*xhat-y{j}) + alpha*norm(P{1}*xhat,1) )
    cvx_end

    x_TV{j} = xhat;

    alpha = 2^(3/2)*(noise_variances{j})/sqrt((sum((Res * x_ls).^2)) / n);

    cvx_begin quiet
        variable xhat(n)
        minimize( norm(F{j}*xhat-y{j}) + alpha*norm(Res*xhat,1) )
    cvx_end

    x_R{j} = xhat;


    nexttile; hold on;
    plot(data_points, x{j},'k--', 'LineWidth',1.5); 
    plot(data_points, y{j}, 'k:', 'LineWidth', 1.5);
    plot(data_points, x_TV{j}, 'b', 'LineWidth',2);
    plot(data_points, x_R{j}, 'color', "#77AC30", 'LineWidth',2);
  
    if j == 1 
        title("Reconstructions");
        xlabel("x (domain)");
        ylabel("Signal value");
        legend("$f_{true}$", ...
            "$y$" ,...
            "${x}_{loc}$",  ...
            "${x}_{R}$", 'Interpreter','latex')
    end
    xlim([-pi, pi]);
    set(gca,'FontSize', 16);

    nexttile; hold on;
    if j == 1 
        title("Pointwise error");
        xlabel("x (domain)");
        ylabel("log_{10}(|error|)");
    end
    plot(data_points, log10(abs(x_TV{j}-x{j})), 'LineWidth',2);
    plot(data_points, log10(abs(x_R{j}-x{j})), 'x', 'LineWidth',2);
    xlim([-pi, pi]);
    set(gca,'FontSize', 16);


    nexttile; hold on;
    plot(data_points, (Res*x_TV{j}), 'LineWidth',2);
    plot(data_points, (Res*x_R{j}), 'x', 'LineWidth',2);
    if j == 1 
        legend('$$R\textit{\textbf{x}}_{TV}$$', '$$R\textit{\textbf{x}}_R$$','interpreter', 'latex'); 
        title("Residual Responses");
        xlabel("x (domain)");
        ylabel("Residual value");
    end
    xlim([-pi, pi]);
    set(gca,'FontSize', 16);
end


%% All visual stuff, plots solutions calculated above
showLegend = false;

if exist('SNR', 'var') && SNR == 20
    showLegend = true;
elseif exist('gamma', 'var') && gamma == 1e-2
    showLegend = true;
elseif exist('r', 'var') && r == 0.1
    showLegend = true;
end

showLegend = false;

j=1;

figure; hold on;
plot(data_points, x{j},'k--', 'LineWidth', 3); 
plot(data_points, y{j}, 'k:', 'LineWidth', 2);
if showLegend
    legend("$$\textit{\textbf{f}}$$", ...
        "$$\textit{\textbf{y}}$$" , 'Interpreter','latex')
end
xlim([-pi, pi]);
set(gca,'FontSize', 20);

figure; hold on;
plot(data_points, x{j},'k--', 'LineWidth', 2); 
plot(data_points, y{j}, 'k:', 'LineWidth', 2);
plot(data_points, x_TV{j}, 'b-x', 'LineWidth',2);
plot(data_points, x_R{j}, 'r-o', 'LineWidth',2);
if showLegend
    legend("$$\textit{\textbf{f}}$$", ...
        "$$\textit{\textbf{y}}$$" ,...
        "$$\textit{\textbf{x}}_{loc}$$",  ...
        "$$\textit{\textbf{x}}_{R}$$", 'Interpreter','latex')
end
xlim([-pi, pi]);
set(gca,'FontSize', 20);

figure; hold on;
plot(data_points, log10(abs(x_TV{j}-x{j})), 'b-x','LineWidth',2);
plot(data_points, log10(abs(x_R{j}-x{j})), 'r-o', 'LineWidth',2);
if showLegend
    legend("$$\textit{\textbf{x}}_{loc}$$",  "$$\textit{\textbf{x}}_R$$", 'Interpreter','latex');
end
xlim([-pi, pi]);
set(gca,'FontSize', 25);


figure; hold on;
plot(data_points, Res*x_TV{j}, 'b-x', 'LineWidth',2);
plot(data_points, Res*x_R{j}, 'r-o', 'LineWidth',2);
if showLegend
    legend("$$R\textit{\textbf{x}}_{loc}$$",  "$$R\textit{\textbf{x}}_R$$", 'Interpreter','latex');
end
xlim([-pi, pi]);
set(gca,'FontSize', 25);


%% Compairison of TV, Residual, and Combination
figure;
tlo = tiledlayout(J,3,'TileSpacing','compact','Padding','compact');

for j = 1:J
    cvx_begin quiet
        variable xhat(n)
        minimize( norm(F{j}*xhat-y{j}) )
    cvx_end

    x_ls = xhat;

    alpha = 2^(3/2)*(noise_variances{j})/sqrt(sum((P{1}*x_ls).^2)/n);

    cvx_begin quiet
        variable xhat(n)
        minimize( norm(F{j}*xhat-y{j}) + alpha*norm(P{1}*xhat,1) )
    cvx_end

    x_TV{j} = xhat;

    alpha = 2^(3/2)*(noise_variances{j})/sqrt((sum((Res * x_ls).^2)) / n);

    cvx_begin quiet
        variable xhat(n)
        minimize( norm(F{j}*xhat-y{j}) + alpha*norm(Res*xhat,1) )
    cvx_end

    x_R{j} = xhat;

    thresh = beta * 1/n; %noise_variances{j} alternative
    x_C{j} = combine_grad_v2(x_TV{j}, x_R{j}, thresh);


    nexttile; hold on;
    plot(data_points, x{j},'k--', 'LineWidth',1.5); 
    plot(data_points, y{j}, 'k:', 'LineWidth', 1.5);
    plot(data_points, x_TV{j}, 'b', 'LineWidth',2);
  
    if j == 1 
        title("TV");
        xlabel("x (domain)");
        ylabel("Signal value");
    end
    xlim([-pi, pi]);
    set(gca,'FontSize', 16);

    nexttile; hold on;
    plot(data_points, x{j},'k--', 'LineWidth',1.5); 
    plot(data_points, y{j}, 'k:', 'LineWidth', 1.5);
    plot(data_points, x_R{j}, 'color', "#77AC30", 'LineWidth',2);
  
    if j == 1 
        title("Residual");
        xlabel("x (domain)");
        ylabel("Signal value");
    end
    xlim([-pi, pi]);
    set(gca,'FontSize', 16);

    nexttile; hold on;
    plot(data_points, x{j},'k--', 'LineWidth',1.5); 
    plot(data_points, y{j}, 'k:', 'LineWidth', 1.5);
    plot(data_points, x_C{j}, 'color', "#77AC30", 'LineWidth',2);
  
    if j == 1 
        title("Combination");
        xlabel("x (domain)");
        ylabel("Signal value");
    end
    xlim([-pi, pi]);
    set(gca,'FontSize', 16);
end

%% define functions
% defining the functions 

%% local differencing
% Constructs matrix B such that (B*x)_j = sum_l w_l (x_{j+1+l} - x_{j-l})
% with periodic indexing


function B = LocalDetector(n, p)
B = zeros(n);  % initialize matrix
q0_p = nchoosek(2*p, p);

for l = 0:p
    w = (-1)^l * nchoosek(2*p + 1, p - l);  % weight

    % shift right: x_{j+1+l}
    shift_fwd = mod((1:n) + 1 + l - 1, n) + 1;  % MATLAB-style wraparound
    B_fwd = zeros(n);
    B_fwd(sub2ind([n n], (1:n)', shift_fwd')) = w;

    % shift left: x_{j-l}
    shift_bwd = mod((1:n) - l - 1, n) + 1;
    B_bwd = zeros(n);
    B_bwd(sub2ind([n n], (1:n)', shift_bwd')) = -w;

    % accumulate both
    B = B + B_fwd + B_bwd;
end

B = B/q0_p;

end

%% contruct_F_undersampling 
%  Function that constructs the forward operator corresponding to 
%  undersampling with parameter ratio r

function F = contruct_F_undersampling(n,r)
    % Define the sampling ratio
    sampling_ratio = 1-r;
    M = round(sampling_ratio * n);

    % Step 1: Randomly select M indices from 1 to N
    random_indices = randperm(n, M);

    % Step 2: Sort the selected indices
    sorted_indices = sort(random_indices);

    % Step 3: Create a logical vector indicating observed variables
    observed = false(n, 1);
    observed(sorted_indices) = true;
    values = double(observed); % Convert logical to double (1s and 0s)

    % Construct the sparse diagonal matrix H
    H = spdiags(values, 0, n, n);
    
    F = eye(n)*H;
end



%% construct_F_deconvolution 
%  Function that constructs the forward operator corresponding to 
%  convolution using a Gaussian kernel with parameter gamma 

function F = construct_F_deconvolution( n, gamma )

    kernel = @(t) exp( -t.^2/(2*gamma^2) )/sqrt(2*pi*gamma^2); % kernel 
    grid = linspace(0, 1, n); % equidistant grid points
    F = kernel( grid-grid' )/n;
    if gamma <5e-3
        F = F/(max(F(:)));
    end
end



%% contruct the concentration factor
% types: trig, exp, loc, polynomial; the filtered version of loc and exp
function sigma = getConcFac(type, p)
    if strcmp(type, 'trig')
        sigma = @(eta) sin(pi*eta)/sinint(pi);

    elseif strcmp(type, 'exp')
        sigma_c = @(eta) exp(1 ./ (p * eta .* (eta - 1)));
        c = 1/quad(@(eta) sigma_c(eta), 0,1);
        sigma = @(eta) c*eta.*exp(1./(p*eta.*(eta-1)));

    elseif strcmp(type, 'loc')
        q0_p = nchoosek(2*p,p);
        sigma = @(eta) (2^(2*p) * eta .* (sin((pi/2)*eta)).^(2*p)) / q0_p;

    elseif strcmp(type, 'poly')
        sigma = @(eta) (2*p+1)*eta.^(2*p+1);
    elseif strcmp(type, 'locFilt')
        q0_p = nchoosek(2*p,p);
        alpha = 32; l=4;
        sigma = @(eta) (exp(-alpha*eta.^l)).*(2^(2*p) * eta .* (sin((pi/2)*eta)).^(2*p)) / q0_p;
    elseif strcmp(type, 'expFilt')
        sigma_c = @(eta) exp(1 ./ (p * eta .* (eta - 1)));
        c = 1/quad(@(eta) sigma_c(eta), 0,1);
        alpha = 32; l= 2;
        sigma = @(eta) c*(exp(-alpha*eta.^l)).*eta.*exp(1./(p*eta.*(eta-1)));

    % Add more elseif blocks for other types as needed

    else
        % --- Handle unknown type ---
        warning('Unknown type: %s. No processing performed.', type);
        return
    end

    nor = quad(@(x)(sigma(x)./x),0,1);
    
    if abs(nor-1)>1e-3
        warning('Wrong normalization: %.4f.', nor);
    end
end




%% construct operator matrix for given concentration factor type
function A = getConMatrix(type,n,p,z)

    if nargin<4
        z = 0;
    end

    % Parameters
    delta_s = 2 * pi / n;
    k = 1 : (n/2);                                % frequency index
    eta = k * delta_s /pi;
    sigma = getConcFac(type,p);
    sin_val = sigma(eta);      
    sin_val(not(isfinite(sin_val))) = 0;

    % Initialize matrix
    A = zeros(n);

    % Loop to fill matrix
    for j = 1:n
        for l = 1:n
            theta1 = k * delta_s * (j - l + 1 - z);
            theta2 = k * delta_s * (j - l     - z);
            kernel = sin_val./k .* (cos(theta1) - cos(theta2));
            A(j, l) = sum(kernel);
        end
    end
end



%% TV_operator
%  Function to compute the one-dimensional TV operator 

function R = TV_operator( n, order )

    e = ones(n,1);
    
    if order == 1
        D = spdiags([e -e], 0:1, n,n);
 	elseif order == 2
       	D = spdiags([-e 2*e -e], 0:2, n,n);
 	elseif order == 3
      	D = spdiags([e -3*e 3*e -e], 0:3, n,n); 
    else 
        error('Desried order not yet implemented!')
    end
    
    R = D(1:n-order,:);
        
end

%% Construct Solution using TV and Residual methods combine
% this function takes two solutions, tv and residual, and combines them
% at each point based on which solution is better. 

function x_com = combine(x_TV, x_R, thresh)
    
    n = length(x_TV);
    x_com = zeros(n, 1);
    abs_diff = abs(x_TV-x_R);

    for i = 1:n
        if abs_diff(i) < thresh
            x_com(i) = x_TV(i);
        else
            x_com(i) = x_R(i);
        end
    end

end

%% Construct Solution using TV and Residual methods combine
% this function takes two solutions, tv and residual, and combines them
% base on the difference in the solution gradient

function x_com = combine_grad(x_TV, x_R, thresh)
    
    n = length(x_TV);
    x_com = zeros(n, 1);
    R_grad = zeros(n, 1);
    TV_grad = zeros(n, 1);
    abs_grad_diff = zeros(n, 1);
    
    for i = 1:n
        if i == 1
            R_grad(i)  = x_R(i+1) - x_R(i);
            TV_grad(i) = x_TV(i+1) - x_TV(i);
        elseif i == n
            R_grad(i)  = x_R(i) - x_R(i-1);
            TV_grad(i) = x_TV(i) - x_TV(i-1);
        else
            R_grad(i)  = (x_R(i+1) - x_R(i-1));   % central diff
            TV_grad(i) = (x_TV(i+1) - x_TV(i-1));
        end
        abs_grad_diff(i) = abs(R_grad(i) - TV_grad(i));
    end

    for i = 1:n
        if abs_grad_diff(i) < thresh
            x_com(i) = x_TV(i);
        else
            x_com(i) = x_R(i);
        end
    end

end

%% Construct Solution using TV and Residual methods combine
% this function takes two solutions, tv and residual, and combines them
% base on the difference in the solution gradient

function x_com = combine_grad_v2(x_TV, x_R, thresh)
    
    n = length(x_TV);
    x_com = zeros(n, 1);
    R_grad = zeros(n, 1);
    TV_grad = zeros(n, 1);
    abs_grad_diff = zeros(n, 1);
    
    for i = 1:n
        if i == 1
            R_grad(i)  = x_R(i+1) - x_R(i);
            TV_grad(i) = x_TV(i+1) - x_TV(i);
        elseif i == n
            R_grad(i)  = x_R(i) - x_R(i-1);
            TV_grad(i) = x_TV(i) - x_TV(i-1);
        else
            R_grad(i)  = (x_R(i+1) - x_R(i-1));   % central diff
            TV_grad(i) = (x_TV(i+1) - x_TV(i-1));
        end
        abs_grad_diff(i) = abs(R_grad(i) - TV_grad(i));
    end

    for i = 1:n
        if abs_grad_diff(i) < thresh
            x_com(i) = (0.8* x_TV(i) + 0.2 * x_R(i));
        else
            x_com(i) = x_R(i);
        end
    end

end
