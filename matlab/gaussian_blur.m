function [Ib] = gaussian_blur(I, wndSize, sigma)
% GAUSSIAN_BLUR Smooth image with symmetric Gaussian filter.
%
%   [Ib] = GAUSSIAN_BLUR(I, wndSize, sigma) produces a filtered image Ib 
%   from I using a square Gaussian kernel with window size wndSize.
%
%   Inputs:
%   -------
%    I        - mxn intensity image.
%    wndSize  - Kernel window size (square, odd number of pixels).
%    sigma    - Standard deviation of Gaussian (pixels, symmetric).
%
%   Outputs:
%   --------
%    Ib  - mxn filtered output image, of same size and class as I.

%--- FILL ME IN ---

% Construct the gaussian smoothing kernel
k = wndSize;
kernel = zeros(k);
for i = 1:k
    for j = 1:k
        kernel(i,j) = exp(-( (i - k/2)^2 + (j - k/2)^2) / (2 * sigma^2));
    end
end
total = sum(sum(kernel));
kernel = kernel / total;

% Convolve the kernel over the image
Ib = conv2(I, kernel, 'same');
if isa(I,'uint8')
    Ib = uint8(round(Ib));
end
if isa(I,'double')
    Ib = double(Ib);
end

%------------------
  
end
