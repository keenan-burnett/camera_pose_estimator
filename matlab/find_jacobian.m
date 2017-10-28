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
%    J  - 2x6 Jacobian matrix (columns are tx, ty, tz, r, p, y).

%--- FILL ME IN ---

R = Ec(1:3,1:3);
t = Ec(1:3,4);
y1 = Wpt - t;
f = K * R' * y1;     % f = distortion free, n = normalized, n = f / z, J = dndP
z = f(3);

[dRdr, dRdp, dRdy] = dcm_jacob_rpy(R);  % Each is 3x3
dRdr = K * dRdr' * y1;
dRdp = K * dRdp' * y1;
dRdy = K * dRdy' * y1;
dfdtheta = [dRdr dRdp dRdy];

dfdt = K * R.' * (-1);

dfdP = [dfdt dfdtheta];     
dzdP = dfdP(3,:);   % P = parameter vector


J = (dfdP / z) - (f * dzdP) / z ^ 2;
J = J(1:2,:);
%------------------

end
