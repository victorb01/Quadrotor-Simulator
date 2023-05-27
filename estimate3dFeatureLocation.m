function [rXIHat,Px] = estimate3dFeatureLocation(M,P)
% estimate3dFeatureLocation : Estimate the 3D coordinates of a feature point
% seen by two or more cameras with known pose.
%
%
% INPUTS
%
% M ---------- Structure with the following elements:
%
% rxArray = 1xN cell array of measured positions of the feature point
% projection on the camera's image plane, in pixels.
% rxArray{i} is the 2x1 vector of coordinates of the feature
% point as measured by the ith camera. To ensure the
% estimation problem is observable, N must satisfy N >= 2 and
% at least two cameras must be non-colinear.
%
% RCIArray = 1xN cell array of I-to-camera-frame attitude matrices.
% RCIArray{i} is the 3x3 attitude matrix corresponding to the
% measurement rxArray{i}.
%
% rcArray = 1xN cell array of camera center positions. rcArray{i} is
% the 3x1 position of the camera center corresponding to the
% measurement rxArray{i}, expressed in the I frame in meters.
%
% P ---------- Structure with the following elements:
%
% sensorParams = Structure containing all relevant parameters for the quad's
% sensors, as defined in sensorParamsScript.m
%
% OUTPUTS
%
%
% rXIHat -------- 3x1 estimated location of the feature point expressed in I
% in meters.
%
% Px ------------ 3x3 error covariance matrix for the estimate rxIHat.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:
%+==============================================================================+

N = length(M.rxArray);
K = P.sensorParams.K;
H = zeros(2*N,4);
R = zeros(2*N);
Rc = P.sensorParams.Rc;
ps = P.sensorParams.pixelSize;

for i = 1:N
    rxVec = ps*M.rxArray{i};
    RCI = M.RCIArray{i};
    rcVec = M.rcArray{i};

    t = -RCI*rcVec;
    P = K*[RCI, t];

    H(2*i-1,:) = (rxVec(1)*P(3,:)-P(1,:));
    H(2*i,:) = (rxVec(2)*P(3,:)-P(2,:));
    R(2*i-1:2*i,2*i-1:2*i) = ps^2*Rc;
end

Hr = H(:,1:3);
z = -H(:,4);
Rinv = inv(R);

Px = inv(Hr'*Rinv*Hr);
rXIHat = Px*Hr'*Rinv*z;

end