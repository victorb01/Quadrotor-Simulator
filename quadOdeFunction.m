function Xdot = quadOdeFunction(t,X,omegaVec,distVec,P)
% quadOdeFunction : Ordinary differential equation function that models
%                   quadrotor dynamics.  For use with one of Matlab's ODE
%                   solvers (e.g., ode45).
%
%
% INPUTS
%
% t ---------- Scalar time input, as required by Matlab's ODE function
%              format.
%
% X ---------- Nx-by-1 quad state, arranged as 
%
%              X = [rI',vI',RBI(1,1),RBI(2,1),...,RBI(2,3),RBI(3,3),omegaB']'
%
%              rI = 3x1 position vector in I in meters
%              vI = 3x1 velocity vector wrt I and in I, in meters/sec
%             RBI = 3x3 attitude matrix from I to B frame
%          omegaB = 3x1 angular rate vector of body wrt I, expressed in B
%                   in rad/sec
%
% omegaVec --- 4x1 vector of rotor angular rates, in rad/sec.  omegaVec(i)
%              is the constant rotor speed setpoint for the ith rotor.
%
%  distVec --- 3x1 vector of constant disturbance forces acting on the quad's
%              center of mass, expressed in Newtons in I.
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
% OUTPUTS
%
% Xdot ------- Nx-by-1 time derivative of the input vector X
%
%+------------------------------------------------------------------------------+

%Unpack inputs
rI = X(1:3);
vI = X(4:6);
RBI = zeros(3,3);
RBI(:) = X(7:15);
omegaB = X(16:18);

%Calculate sum of Forces
sumF = zeros(3,1);
for i = 1:4
    sumF = sumF + P.quadParams.kF(i)*omegaVec(i)^2*[0; 0; 1];
end

%Calculate sum of Torques
sumN = zeros(3,1);
for i = 1:4
    sumN = sumN - P.quadParams.omegaRdir(i)*P.quadParams.kN(i)*omegaVec(i)^2*[0; 0; 1];
end

%Calculate the r cross F components of torque

sumNF = zeros(3,1);
for i = 1:4
    sumNF = sumNF + P.quadParams.kF(i)*omegaVec(i)^2*crossProductEquivalent(P.quadParams.rotor_loc(:,i))*[0; 0; 1];
end

%Calculate the derivatives
rIdot = vI;
vIdot = [0; 0; -P.constants.g] + transpose(RBI)*sumF/P.quadParams.m + distVec/P.quadParams.m;
RBIdot = -crossProductEquivalent(omegaB)*RBI;
RBIdotVec = RBIdot(:);
omegaBdot = inv(P.quadParams.Jq)*(sumN + sumNF - crossProductEquivalent(omegaB)*P.quadParams.Jq*omegaB);

%Make the final output
Xdot = [rIdot; vIdot; RBIdotVec; omegaBdot];
end