function [pt] = saddle_point(I)
% SADDLE_POINT Locate saddle point in an image patch.
%
%   [pt] = SADDLE_POINT(I) finds the subpixel center of a cross-junction in the 
%   image patch I, by blurring the patch, fitting a hyperbolic paraboloid to 
%   it, and then finding the critical point of that paraboloid.
%
%   Note that the location of 'p' is relative to (0.5, 0.5) at the upper left 
%   corner of the patch, i.e., the pixels are treated as covering an area of
%   one unit square.
%
%   Inputs:
%   -------
%    I  - mxn image patch (grayscale, double or integer class).
%
%   Outputs:
%   --------
%    pt  - 2x1 subpixel location of saddle point in I (x, y coords).
%
% References:
%
%   L. Lucchese and S. K. Mitra, "Using Saddle Points for Subpixel Feature
%   Detection in Camera Calibration Targets," in Proc. Asia-Pacific Conf.
%   Circuits and Systems (APCCAS'02), vol. 2, (Singapore), pp. 191-195,
%   Dec. 2002.

%--- FILL ME IN ---

[m,n] = size(I);
sampleRate = 100;
I = im2double(I);

% Do a bicubic interpolation over the image patch
[X,Y] = meshgrid(1:n,1:m);
[Xfine,Yfine] = meshgrid(1:1/sampleRate:n,1:1/sampleRate:m);
Ifine = interp2(X,Y,I,Xfine,Yfine, 'cubic');
[M,N] = size(Ifine);

% Get the gradients of the image
kx = [-1 0 1; -1 0 1; -1 0 1];
ky = kx';
dfdy = conv2(Ifine,ky,'same');
dfdx = conv2(Ifine,kx,'same');
d2fdy2 = conv2(dfdy,ky,'same');
d2fdx2 = conv2(dfdx,kx,'same');
d2fdxdy = conv2(dfdy,kx,'same');

% The Hessian (H) and Magnitude (Mag) are used to determine saddle points 
H = d2fdx2 .* d2fdy2 - d2fdxdy;
H = H(3:end-3,3:end-3);
Mag = sqrt(dfdx.^2 + dfdy.^2);
Mag = Mag(3:end-3,3:end-3);

% Threshold on Magnitude and the Hessian value
mag_t = max(Mag(:)) / 16;
saddles = Mag < mag_t & H < 10;

% Find 'center-most' point among points produced
[row,col] = find(saddles);
r = (row - M/2).^2 + (col - N/2).^2;
[~,index] = min(r);
y = row(index) + 3;
x = col(index) + 3;
pt = [x;y] / sampleRate;
%------------------
  
end
