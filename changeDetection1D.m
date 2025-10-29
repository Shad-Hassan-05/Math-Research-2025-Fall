
% Free parameters of the problem 
n = 128; % number of (equidistant) data points on [0,1] 
SNR = 15; % SNR used to compute variance of the iid Gaussian noise added to the measurements
% gamma = 1e-2

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
f = f1(jumpLoc);

% creates noisy function y
x = f(data_points); % function values at grid points 
F = eye(n);
y = F*x;
noise_variances = norm(x)^2/n/(10^(SNR/10));
noise = sqrt(noise_variances)*randn(size(F,1),1);
y = y + noise; % noisy data 
snr = 10*log10(norm(x)^2/n/noise_variances);

% plot noisy y
figure;
plot(data_points, x, 'LineWidth',2); hold on;
plot(data_points, y,'LineWidth',1.5); hold on;
legend("${f}$", "${y}$", 'fontsize', 16, 'Interpreter', 'latex');
title("True vs Noisy Signal");
xlabel("x (domain)");
ylabel("Signal value"); 

%% Getting two reconstruction images using different methods
% Grid and frequency setup
% builds the residual operator. S is gloabal fourier, T is the local diff.
% this is what we compair to the TV method. 
p = 0;

S1 = getConMatrix('loc',n,p + 1);
S2 = getConMatrix('exp',n,p + 3);

edgeMap1 = S1*y+S1'*y;
edgeMap2 = S2*y+S2'*y;

figure;
plot(data_points, edgeMap1); hold on;
plot(data_points, edgeMap2); hold on;
legend("${S1}$", "${S2}$", 'fontsize', 16, 'Interpreter', 'latex');
title("edgeMap1 vs edgeMap2");

%% Use 1D CDD to differenciate change
% using code from ChanfeDectection_Old_Ash file
m = n/2;
k1 = (m-1)/2;
k2 = 3*(m-1)/2;

% neighborhood size
nbr_size = 10;

% Store test statistic here
gamma_CCD = zeros(n,1); 

y1 = edgeMap1;
y2 = edgeMap2;

% Loop which tells us at which indices to compute the test statistic
for idx_k = nbr_size:(n-1)-nbr_size

    % Extract appropriate neighborhood of test signals
    y1_nbr = y1((idx_k+1)-nbr_size:(idx_k+1)+nbr_size);      
    y2_nbr = y2((idx_k+1)-nbr_size:(idx_k+1)+nbr_size);


    %  Compute test statistic
    gamma_CCD(idx_k+1) = 1 - 2*abs(y2_nbr'*y1_nbr)/(sum( abs(y1_nbr).^2 + abs(y2_nbr).^2 ));
    
end

figure;
plot(data_points, gamma_CCD, 'LineWidth', 1.8);
xlabel('x'); 
ylabel('\gamma_{CCD}');
title('1D CCD Test Statistic');
grid on;