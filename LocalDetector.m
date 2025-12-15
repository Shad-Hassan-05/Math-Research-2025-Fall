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