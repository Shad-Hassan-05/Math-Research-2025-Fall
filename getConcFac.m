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

function [outputArg1,outputArg2] = untitled2(inputArg1,inputArg2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end