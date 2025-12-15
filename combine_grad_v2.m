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