function [xkp1] = f_dynamics(xk,uk,vk,delt,RBIHatk,P)
% f_dynamics : Discrete-time dynamics model for quadcopter.
%
% INPUTS
%
% xk --------- 15x1 state vector at time tk, defined as 
% 
%              xk = [rI', vI', e', ba', bg']'
%
%              where all corresponding quantities are identical to those
%              defined for E.statek in stateEstimatorUKF.m and where e is the
%              3x1 error Euler angle vector defined such that for an estimate
%              RBIHat of the attitude, the true attitude is RBI =
%              dRBI(e)*RBIHat, where dRBI(e) is the DCM formed from the error
%              Euler angle vector e.
%
% uk --------- 6x1 IMU measurement input vector at time tk, defined as
%
%              uk = [omegaBtilde', fBtilde']'
%
%              where all corresponding quantities are identical to those
%              defined for M in stateEstimatorUKF.m.
%
% vk --------- 12x1 process noise vector at time tk, defined as
%
%              vk = [vg', vg2', va', va2']'
%
%              where vg, vg2, va, and va2 are all 3x1 mutually-independent
%              samples from discrete-time zero-mean Gaussian noise processes.
%              These represent, respectively, the gyro white noise (rad/s),
%              the gyro bias driving noise (rad/s), the accelerometer white
%              noise (m/s^2), and the accelerometer bias driving noise
%              (m/s^2).
%
% delt ------- Propagation interval, in seconds.
%
% RBIHatk ---- 3x3 attitude matrix estimate at time tk.
%
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
%  sensorParams = Structure containing sensor parameters, as defined in
%                 sensorParamsScript.m
%
%
% OUTPUTS
%
% xkp1 ------- 15x1 state vector propagated to time tkp1
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: 
%+==============================================================================+  

if(abs(delt - P.sensorParams.IMUdelt) > 1e-9)
  error('Propagation time must be same as IMU measurement time');
end

%Initialise given variables
rIk = xk(1:3);
vIk = xk(4:6);
ek = xk(7:9);
bak = xk(10:12);
bgk = xk(13:15);
omegaBtildek = uk(1:3);
fBtildek = uk(4:6);
vgk = vk(1:3);
vg2k = vk(4:6);
vak = vk(7:9);
va2k = vk(10:12);
RBIk = euler2dcm(ek)*RBIHatk;
phi = ek(1);
theta = ek(2);
g = P.constants.g;
e3 = [0 0 1]';
alphaa = P.sensorParams.alphaa;
alphag = P.sensorParams.alphag;

aIk = (RBIk')*(fBtildek - bak - vak) - g*e3;
rIkp1 = rIk + delt*vIk + 0.5*delt^2*aIk;
vIkp1 = vIk + delt*aIk;
gamma = [cos(theta), 0, sin(theta);
    sin(phi)*sin(theta)/cos(phi), 1, -cos(theta)*sin(phi)/cos(phi);
    -sin(theta)/cos(phi), 0, cos(theta)/cos(phi)];
ekdot = gamma*(omegaBtildek - bgk - vgk);
ekp1 = ek + delt*ekdot;
bakp1 = alphaa*bak + va2k;
bgkp1 = alphag*bgk + vg2k;

xkp1 = [rIkp1; vIkp1; ekp1; bakp1; bgkp1];

end