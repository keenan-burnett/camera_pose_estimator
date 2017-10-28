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

threshold = 10;
E = Eg;
R = E(1:3,1:3);
t = E(1:3,4);
q = rpy_from_dcm(R);
P = [t ; q];
K = [564.9 0 337.3; 0 564.3 226.5; 0 0 1];
damp = eye(6) * 0.01;
cost = 0;
lastCost = 1e50;

while (lastCost - cost) > threshold
    % Get the rotation and translation from homogeneous pose
    R = E(1:3,1:3);
    t = E(1:3,4);
    z = t(3);
    
    % Calculate A and B
    A = 0;
    b = 0;
    for i = 1:48
        xi = Wpts(:,i);
        J_xi = find_jacobian(K, E, xi);
        A = A + J_xi' * J_xi;
        r_i = Ipts(:,i) - (K * R' * (xi - t) / z);
        b = b + J_xi' * r_i;
    end
    
    % Update P
    deltaP = b / (A + damp);
    P = P + deltaP;
    
    lastCost = cost;
    cost = 0;
    % Calculate cost
    for i = 1:48
        xi = Wpts(:,i);
        J_xi = find_jacobian(K, E, xi);
        r_i = Ipts(:,i) - (K * R' * (xi - t) / z);
        cost = cost + norm(J_xi * deltaP - r_i) ^ 2;
    end
        
    % Update pose matrix
    E(1:3,1:3) = dcm_from_rpy(P(4:6));
    E(1:3,4) = P(1:3);
    
end
    
    
        
        
        
    
end

%------------------

end
