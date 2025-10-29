clear, close all;

% to do: Use Tv : w = Tx + TT' and compair w and shifted w in CCD

% Define input parameters for toy image generation
n = 200;
J = 3;
radius1 = pi*.3;
radius2 = pi*.7;

% construct functions
[pointx, pointy, f] = generate_images(J, n, radius1, radius2);

% get image made by function
im1 = f{3}; %+ 0.2 * randn(size(f{3}));
im2 = f{3};% + 0.5 * randn(size(f{3}));

dy = 1;
dx = 1;

% tv image and shifted tv image
tv_im = tv(im1);
shifted_tv = circshift(tv_im, [dy, dx]);

% shifted image
shifted_im = circshift(im1, [dy, dx]);

% con matrix
s1 = getConMatrix('poly', n, 2);
s2 = getConMatrix('loc', n, 0);

edgeMap1 = s1*im1+im1*s1';
edgeMap2 = s2*im1+im1*s2';

figure; imagesc(im1);

figure;
imagesc(edgeMap1);

figure; 
imagesc(edgeMap2);

idx = n/2;
figure;
hold on;
% plot(im1(idx,:));
plot(edgeMap1(idx,:));
plot(edgeMap2(idx,:));


% define parameters for CCD function
CCDwinSize = 10;

%% Process image and compute change map in each pixel using CCDbasic.m

% Perform CCD using the MLE method
CCDimage = CCDbasic(edgeMap1,edgeMap2,CCDwinSize,'MLE',false);
% the last argument determines whether to implement low-RCS masking

% Display the CCD results
figure(1);
imagesc(1 - CCDimage,[0 1]);
axis image
title('Coherent Change Detection (MLE)');
colormap gray
colorbar

%% TV function for 2D

function tv = tv(u)
    u = double(u);
    [ux, uy] = gradient(u);
    tv = sqrt(ux.^2 + uy.^2);
end


%% function to create the toy images
% taken from script_imag_gen.m 

function [pointx, pointy, f] = generate_images(J, n, r1, r2)
% GENERATE_IMAGE_MATRICES creates three image matrices with common edges.
%
%
% Syntax:
%   [x, y, f] = generate_images(n, r1, r2, J)
%
% Inputs:
%   n  - Integer, specifies the image dimension (n x n).
%   r1 - Double, the radius of the inner circular edge.
%   r2 - Double, the radius of the outer circular edge.
%   J  - Integer, the number of images to generate (up to 3).
%
% Outputs:
%   x  - n x 1 vector of x-coordinates [0, 1).
%   y  - n x 1 vector of y-coordinates, flipped version of x.
%   f  - 1xJ cell array containing the generated image matrices.
%

% 1. Sampling and coordinate setup for the [-pi, pi] domain
xs = -pi; xe = pi;
xt = linspace(xs, xe, n+1)';
pointx = xt(1:end-1);
pointy = flipud(pointx);
[xx, yy] = meshgrid(pointx, pointy);

% 2. Define the common edges using radii
radius_map = sqrt(xx.^2 + yy.^2);

% Create boolean masks for each of the three regions
inner_mask = radius_map < r1;
middle_mask = (radius_map >= r1) & (radius_map < r2);
outer_mask = radius_map >= r2;

% 3. Generate matrices and assign to output cell array
f = cell(1, J);

% --- Generate f1 ---
f1 = zeros(n);
periodic_poly = (cos(xx) + 1).^2 + (cos(yy) + 1).^2;
f1(inner_mask) = 3 + periodic_poly(inner_mask);
f1(outer_mask) = sin(xx(outer_mask));
f1(middle_mask) = 4 + sin(4*xx(middle_mask)) + sin(4*yy(middle_mask));
f{1} = f1;

% --- Generate f2 if J >= 2 ---
if J >= 2
    f2 = zeros(n);
    f2(inner_mask) = sin(6 * xx(inner_mask));
    f2(middle_mask) = -0.3 * sin(6 * yy(middle_mask));
    f2(outer_mask) = sin(-xx(outer_mask) + pi);
    f{2} = f2;
end

% --- Generate f3 if J >= 3 ---
if J >= 3
    f3 = zeros(n);
    % This expression is correct because the negative sign from the transformation
    % `cos(xx+pi) = -cos(xx)` is squared, becoming positive.
    periodic_radius_map = sqrt((2*pi*cos(xx)).^2 + (2*pi*cos(yy)).^2);
    
    f3(inner_mask) = cos(2 * radius_map(inner_mask));
    f3(middle_mask) = cos(4 * radius_map(middle_mask));

    % This term was corrected by adding the negative sign from the cos(x+pi) identity.
    f3(outer_mask) = 5; % * cos(xx(outer_mask)); % <-- CORRECTED
    
    f{3} = f3;
end
end