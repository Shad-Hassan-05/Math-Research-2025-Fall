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