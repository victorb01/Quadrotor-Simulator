function [eak] = voltageConverter(Fk,NBk,P)
% voltageConverter : Generates output voltages appropriate for desired
%                    torque and thrust.
%
%
% INPUTS
%
% Fk --------- Commanded total thrust at time tk, in Newtons.
%
% NBk -------- Commanded 3x1 torque expressed in the body frame at time tk, in
%              N-m.
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
% eak -------- Commanded 4x1 voltage vector to be applied at time tk, in
%              volts. eak(i) is the voltage for the ith motor.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+

%Set up alpha and beta
alpha = 1;
beta = 0.9;

%Calculate k_T
k_T = P.quadParams.kN./P.quadParams.kF;

%Create G matrix
G = [ones(1,4); P.quadParams.rotor_loc(2,:); -P.quadParams.rotor_loc(1,:); -k_T'.*P.quadParams.omegaRdir];

%Calculate Fmax
omegaMax = min(P.quadParams.cm)*P.quadParams.eamax;
Fmax = min(P.quadParams.kF)*omegaMax^2;

Ginv = inv(G);

%Find vector of rotor forces
F = Ginv*[min(Fk, 4*beta*Fmax); alpha*NBk];
while(max(F) > Fmax)
    alpha = alpha*0.99;
    F = Ginv*[min(Fk, (4*beta*Fmax)); alpha*NBk];
end
F(find(F < 0)) = 0;

%Find voltage vector
eak = sqrt(F./P.quadParams.kF)./P.quadParams.cm;

end