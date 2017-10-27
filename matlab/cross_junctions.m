function [Ipts] = cross_junctions(I, boundPoly, Wpts)
% CROSS_JUNCTIONS Find cross-junctions in image with subpixel accuracy.
%
%   [Ipts] = CROSS_JUNCTION(I, boundPoly, Wpts) locates a series of cross- 
%   junction points on a planar calibration target, where the target is
%   bounded in the image by the specified 4-sided polygon. The number of
%   cross-junctions identified should be equal to the number of world points.
%
%   Note also that the world and image points must be in *correspondence*,
%   that is, the first world point should map to the first image point, etc.
%
%   Inputs:
%   -------
%    I          - Image (grayscale, double or integer class).
%    boundPoly  - 2x4 array defining bounding polygon for target (clockwise).
%    Wpts       - 3xn array of world points (in 3D, on calibration target).
%
%   Outputs:
%   --------
%    Ipts  - 2xn array of cross-junctions (x, y), relative to the upper left
%            corner of I.

%--- FILL ME IN ---

I = double(I);

% Create an polygon mask around the checkerboard
bw = mask(I, boundPoly);

% Perform Harris Corner Detection
kx = [-1 0 1; -2 0 2; -1 0 1];
ky = kx';
Iy = conv2(I,ky,'same');
Ix = conv2(I,kx,'same');
Ix2 = gaussian_blur(Ix.^2,7,2);
Iy2 = gaussian_blur(Iy.^2,7,2);
Ixy = gaussian_blur(Ix.*Iy,7,2);
[h,w] = size(I);
C = zeros(h,w);
f = find(bw);
[bwRows,bwCols] = ind2sub([h,w],f);
[len,~] = size(bwRows);
for i = 1:len
    u = bwCols(i);
    v = bwRows(i);
    M = [Ix2(v,u) Ixy(v,u); Ixy(v,u) Iy2(v,u)];
    C(v,u) = det(M) - 0.04 * (trace(M))^2;
end

% Threshold on harris corner value
threshold = 0.1;
corners = C > threshold*max(C(:));
[row,~] = find(corners);
[len,~] = size(row);
while len > 2500
    threshold = threshold + 0.01;
    corners = C > threshold*max(C(:));
    [row,~] = find(corners);
    [len,~] = size(row);
end

% Erode the points resulting from the threshold
kernel = ones(5);
erode = conv2(corners,kernel,'same');
maximum = max(erode(:));
corners = erode >= maximum;

% Reduce clusters to single points
[row,col] = find(corners);
[len,~] = size(row);
i = 1;
j = 1;
cRow(j,1) = row(i);
cCol(j,1) = col(i);
j = j + 1;
for i = 1:len
    [R,~] = findNearestNeighborDistance(cRow,cCol,row(i),col(i));
    if R > 5
        cRow(j,1) = row(i);
        cCol(j,1) = col(i);
        j = j + 1;
    end
end

% Perform saddle point detection
corners = [cCol cRow]';
for i = 1:48
    x = corners(1,i);
    y = corners(2,i);
    snip = I(y-2:y+4,x-2:x+4);
    pt = saddle_point(snip);
    corners(1,i) = pt(1) + x - 2;
    corners(2,i) = pt(2) + y - 2;
end

% Sort the points in row-major order
Ipts = sort_corners(corners);

% Show: x junctions after subpixel estimation
imshow(uint8(I));
hold on
plot(corners(1,:)', corners(2,:)', 'r+');
hold off
  
%------------------
end

% Returns a mask delimited by the boundpoly given. 
% Assumes 4 points clockwise from top-left.
function bw = mask(Image, boundpoly)
    [h,w] = size(Image);
    x2 = max(boundpoly(1,:));
    x1 = min(boundpoly(1,:));
    y1 = min(boundpoly(2,:));
    y2 = max(boundpoly(2,:));
    X = (x1 + x2)/2;
    X = round(X);
    Y = (y1 + y2)/2;
    Y = round(Y);
    pts = boundpoly;
    pts(1,:) = pts(1,:) - X;
    pts(2,:) = pts(2,:) - Y;
    % Shrink ROI to exclue the outsides of the checkerboard
    pts = pts * 0.85;
    pts(1,:) = pts(1,:) + X;
    pts(2,:) = pts(2,:) + Y;

    bw = ones(h,w) * 255;
    M1 = [pts(2,1) pts(2,2)] / [pts(1,1) pts(1,2); 1 1];     % Top line  y = mx + b
    M2 = [pts(1,2) pts(1,3)] / [pts(2,2) pts(2,3); 1 1];    % Right line x = my + b
    M3 = [pts(2,3) pts(2,4)] / [pts(1,3) pts(1,4); 1 1];   % Bottom line y = mx + b
    M4 = [pts(1,4) pts(1,1)] / [pts(2,4) pts(2,1); 1 1];    % Left line x = my + b

    % Mask out the points outside the polygon created by the four lines above
    for i = 1:h
        for j = 1:w
            if (i + 2 < M1 * [j; 1]) || (i > M3 * [j; 1]) || (j-1 > M2 * [i; 1]) ...
                    || (j < M4 * [i; 1])
                bw(i,j) = 0;
            end
        end
    end

end

% finds the nearest neighbor in (rows,cols) to (y,x)
function [R,index] = findNearestNeighborDistance(rows,cols,y,x)
    [len,~] = size(rows);
    R = inf;
    index = 0;
    for i = 1:len
        r = (rows(i) - y)^2 + (cols(i) - x)^2;
        if r < R
            R = r;
            index = i;
        end
    end
    R = sqrt(R);
end

% Sorts 'corners' into row-major order
function sorted = sort_corners(corners)
    sorted = zeros(2,48);
    for i = 1:6
        rows = corners(2,:)';
        cols = corners(1,:)';
        r = rows.^2 + cols.^2;
        [~,index] = min(r);
        p = zeros(2,8);
        
        
        p(:,1) = corners(:,index);
        corners(:,index) = [];
        
        for j = 2:8
            p_0 = p(:,j-1);
            rows = corners(2,:)';
            cols = corners(1,:)';
            [~,index1] = findNearestNeighborDistance(rows,cols,p_0(2,1),p_0(1,1));
            p_1 = corners(:,index1);
            corners2 = corners;
            corners2(:,index1) = [];
            rows2 = corners2(2,:)';
            cols2 = corners2(1,:)';
            [~,index2] = findNearestNeighborDistance(rows2,cols2,p_0(2,1),p_0(1,1));
            if i < 6
                p_2 = corners2(:,index2);
                if p_2(1) > p_1(1)
                    p(:,j) = p_2;
                    [~,index2] = findNearestNeighborDistance(rows,cols,p_2(2,1),p_2(1,1));
                    corners(:,index2) = [];
                else
                    p(:,j) = p_1;
                    corners(:,index1) = [];
                end
            else
                p(:,j) = p_1;
                corners(:,index1) = [];
            end
        end
        sorted(:,((i - 1) * 8 + 1): i * 8) = p;    
    end
end
    