clear, close all

% Free parameters of the problem 
n = 256; % number of (equidistant) data points on [0,1] 
SNR = 10; % SNR used to compute variance of the iid Gaussian noise added to the measurements
% gamma = 1e-2

%% Set up the model 

 %(v(1) <= x & x< v(2)) .*( -6/pi ) +...
f1 = @(v) @(x) (x <v(1)) .* (x + pi+1) +...
    (v(1) <= x & x< v(2)) .*(-1/2*sin(6*x)) +...
    (x >= v(2)) .*(1); %sin(-x+pi), -1/2*sin(6*2)

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

%% Plotting noisy signal 

% plot noisy y
figure;
plot(data_points, x, 'LineWidth',2); hold on;
plot(data_points, y,'LineWidth',1.5); hold on;
legend("${f}$", "${y}$", 'fontsize', 16, 'Interpreter', 'latex');
title("True vs Noisy Signal");
xlabel("x (domain)");
ylabel("Signal value"); 

%% get residual solutions.
% changing alpha values

cvx_begin quiet
   variable xhat(n)
   minimize( norm(F*xhat-y) )
cvx_end

x_ls = xhat;

p = 0;

r = 9; % 2^(3/2)*(noise_variances)/sqrt((sum((Res * x_ls).^2)) / n)
%alphas = linspace(0.2, 2.4, r);
alphas = logspace(-.6, .6, r);

residual_norm = zeros(1, r);
reg_norm      = zeros(1, r);

R = zeros(n, r);
Res = zeros(n, r);

for k=1:r

    T = LocalDetector(n,p);
    S = getConMatrix('loc',n,p,1/4);

    T_1 = LocalDetector(n, p+2);

    Rs = (T-S);

    alpha = alphas(k);  %* 2^(3/2)*(noise_variances)/sqrt((sum((Rs * x_ls).^2)) / n);

    cvx_begin quiet
        variable xhat(n)
        minimize( norm(F*xhat-y) + alpha*norm(Rs*xhat,1) )
    cvx_end

    Res(:, k) = xhat;

    alpha = alphas(k) * 2^(3/2)*(noise_variances)/sqrt((sum((T_1 * x_ls).^2)) / n);

    cvx_begin quiet
        variable xhat(n)
        minimize( norm(F*xhat-y) + alpha*norm(T_1*xhat,1) )
    cvx_end

    R(:, k) = xhat;

end

figure; tiledlayout(3, ceil(r/3));  % 3 rows of tiles
for k = 1:r
    nexttile;
    plot(data_points, x, 'k--', 'LineWidth',2); hold on;
    plot(data_points, Res(:,k),'LineWidth',1.5);
    title(sprintf('\\alpha = %.3g', alphas(k)));
    grid on; axis tight;
end

figure
plot(data_points,Res)

% figure;
% loglog(residual_norm, reg_norm, 'o-','LineWidth',1.2);
% xlabel('Residual norm ||F x - y||');
% ylabel('Regularization norm ||Res x||_1');
% title('L-curve');
% grid on;

q=1;

TV = zeros(n, r);
%alphas_tv = linspace(0.2, 2.4, r);
alphas_tv = logspace(-.6, .6, r);

for k=1:r

    P = TV_operator( n, q);

    alpha = alphas_tv(k) * 2^(3/2)*(noise_variances)/sqrt((sum((P * x_ls).^2)) / n);;

    cvx_begin quiet
        variable xhat(n)
        minimize( norm(F*xhat-y) + alpha*norm(P*xhat,1) )
    cvx_end

    TV(:, k) = xhat;

end

figure; tiledlayout(3, ceil(r/3));  % 3 rows of tiles
for k = 1:r
    nexttile;
    plot(data_points, x, 'k--', 'LineWidth',2); hold on;
    plot(data_points, TV(:,k),'LineWidth',1.5);
    title(sprintf('\\alpha = %.3g', alphas(k)));
    grid on; axis tight;
end

%% calvulate the Var for each data point between the different alphas
% var = zeros(size(y));
% 
% for i = 1:n
%     for k = 1:r-1
%         var(i) = var(i) + norm(TV(i,k) - TV(i,k+1));
%     end
% end

VarTV = var(TV,[],2);

figure; hold on;
plot(data_points, VarTV, 'LineWidth',2); hold on;
legend("${Var}$", 'fontsize', 16, 'Interpreter', 'latex');
 axis tight;
title("Var Plot");
% ylim([0,0.2])
% xlabel("x (domain)");
% ylabel("abs diff value"); 
% 

VarR = var(Res,[],2);
Err = mean(abs(Res-repmat(x, 1, r)), 2);

% figure;
plot(data_points, VarR, 'LineWidth',2); hold on;
plot(data_points, Err, 'k--', 'LineWidth',2); hold on;
legend("${Var}$", 'fontsize', 16, 'Interpreter', 'latex');
 axis tight;
 % ylim([0,0.2])
title("Var Plot");
xlabel("x (domain)");
ylabel("abs diff value"); 

%% detecting constant sections using point wise difference.
sol = zeros(size(y));
thresh = noise_variances;
large_alph = Res(:,5);
small_alph = Res(:,3); 
diff = VarTV - VarR;

for i=1:n

    if diff(i) < thresh
        if diff(i) > 0 
            sol(i) = large_alph(i);
        else
            sol(i) = small_alph(i);
        end
    else
        sol(i) = small_alph(i);
    end
end

figure; hold on;
plot(data_points, sol, 'LineWidth',2); hold on;
% plot(data_points, y,'LineWidth',1.5); hold on;
plot(data_points, x,'LineWidth',1.5); hold on;
legend("${sol}$", 'fontsize', 16, 'Interpreter', 'latex');
 axis tight;
title("Mixed Solution");

%% over a interval 
% 
% 
% p = 3;
% T = LocalDetector(n,p);
% 
% q=1;
% P = TV_operator( n, q);  
% 
% 
% alpha = 2^(3/2)*(noise_variances)/sqrt(sum((P*x_ls).^2)/n);
% 
% cvx_begin quiet
%   variable xhat(n)
%   minimize(norm(F*xhat-y) + alpha*norm(P*xhat,1))
% cvx_end
% 
% x_TV = xhat;
% 
% TV_edges = T*x_TV;
% TV_edge_abs = abs(TV_edges);
% 
% figure;
% plot(TV_edges, 'LineWidth',2); hold on;
% plot(TV_edge_abs, 'LineWidth',1.5); hold on;
% legend("${TV_edges}$","${TV_edges_abs}$", 'fontsize', 16, 'Interpreter', 'latex');
% title("TV edge map and Abs TV edge map");
% 
% thresh = 1;
% 
% jumps = [];   % dynamic list okay here
% 
% for j = 1:n
%     if TV_edge_abs > thresh
%         jumps(end+1) = j;   % append index
%     end
% end
% 
% % Now we know how many jumps:
% intervals = zeros(length(jumps)-1,2);
% 
% for k = 1:length(jumps)-1
%     bounds = [1, jumps, n+1];               % note n+1 so last becomes [jm, n]
%     intervals = [bounds(1:end-1)', bounds(2:end)'-1];
% end

%% Add zero padding to the function on both ends at least 5 on each side. 

