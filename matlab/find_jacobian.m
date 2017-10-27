function [J] = find_jacobian(K, Ec, Wpt)
%  FIND_JACOBIAN Determine Jacobian for NLS camera pose optimization.
%
%   [J] = FIND_JACOBIAN(K, Ec, Wpt) computes the Jacobian of image plane point
%   with respect to the current camera pose estimate and given a single world
%   point.
%
%   Inputs:
%   -------
%    K    - 3x3 camera intrinsic calibration matrix.
%    Ec   - 4x4 homogenous pose matrix, current guess for camera pose.
%    Wpt  - 3x1 world point on calibration target (one of n).
%
%   Outputs:
%   --------
%    J  - 2x6 Jacobian matrix (columns are tx, ty, tz, r, p, q).

%--- FILL ME IN ---

% Code goes here...

%------------------

end
