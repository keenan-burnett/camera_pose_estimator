function [E] = pose_estimate_nlopt(Eg, Ipts, Wpts)
%  POSE_ESTIMATE_NLOPT Estimate camera pose from 2D-3D correspondences via NLS.
%
%   [E] = POSE_ESTIMATE_NLOPT(Eg, Ipts, Wpts) performs a nonlinear least squares 
%   optimization procedure to determine the best estimate of the camera pose in 
%   the calibration target frame, given 2D-3D point correspondences.
%
%   Inputs:
%   -------
%    Eg    - 4x4 homogenous pose matrix, initial guess for camera pose.
%    Ipts  - 2xn array of cross-junction points (with subpixel accuracy).
%    Wpts  - 3xn array of world points (one-to-one correspondence with image).
%
%   Outputs:
%   --------
%    E  - 4x4 homogenous pose matrix, estimate of camera pose in target frame.

%--- FILL ME IN ---

% Code goes here...

%------------------

end
