function NBk = attitudeController(R,S,P)
% attitudeController : Controls quadcopter toward a reference attitude
%
%
% INPUTS
%
% R ---------- Structure with the following elements:
%
%       zIstark = 3x1 desired body z-axis direction at time tk, expressed as a
%                 unit vector in the I frame.
%
%       xIstark = 3x1 desired body x-axis direction, expressed as a
%                 unit vector in the I frame.
%
% S ---------- Structure with the following elements:
%
%        statek = State of the quad at tk, expressed as a structure with the
%                 following elements:
%                   
%                  rI = 3x1 position in the I frame, in meters
% 
%                 RBI = 3x3 direction cosine matrix indicating the
%                       attitude
%
%                  vI = 3x1 velocity with respect to the I frame and
%                       expressed in the I frame, in meters per second.
%                 
%              omegaB = 3x1 angular rate vector expressed in the body frame,
%                       in radians per second.
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
%
% OUTPUTS
%
% NBk -------- Commanded 3x1 torque expressed in the body frame at time tk, in
%              N-m.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+

%Create gain matrices
K = diag([1 0.95 0.21]);
K_d = diag([0.25 0.25 0.12]);

%Create desired rotation matrix
b = cross(R.zIstark, R.xIstark);
if b == 0
    b = [0 0 0]';
else
    b = b/norm(b);
end
a = cross(b, R.zIstark);
a = a/norm(a);
RBIstar = [a, b, R.zIstark]';

%Create error direction-cosine matrix
R_e = RBIstar*(S.statek.RBI');

%Create e vector
e_E = [(R_e(2,3)-R_e(3,2)); (R_e(3,1)-R_e(1,3)); (R_e(1,2)-R_e(2,1))];

%Find NBk
NBk = K*e_E - K_d*S.statek.omegaB + crossProductEquivalent(S.statek.omegaB)*(P.quadParams.Jq*S.statek.omegaB);

end