function [rpGtilde,rbGtilde] = gnssMeasSimulator(S,P)
% gnssMeasSimulator : Simulates GNSS measurements for quad.
%
%
% INPUTS
%
% S ---------- Structure with the following elements:
%
%        statek = State of the quad at tk, expressed as a structure with the
%                 following elements:
%                   
%                  rI = 3x1 position of CM in the I frame, in meters
% 
%                 RBI = 3x3 direction cosine matrix indicating the
%                       attitude of B frame wrt I frame
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the
%                 quad's sensors, as defined in sensorParamsScript.m 
%
%
% OUTPUTS
%
% rpGtilde --- 3x1 GNSS-measured position of the quad's primary GNSS antenna,
%              in ECEF coordinates relative to the reference antenna, in
%              meters.
%
% rbGtilde --- 3x1 GNSS-measured position of secondary GNSS antenna, in ECEF
%              coordinates relative to the primary antenna, in meters.
%              rbGtilde is constrained to satisfy norm(rbGtilde) = b, where b
%              is the known baseline distance between the two antennas.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+

%Initialise given variables
rI = S.statek.rI;
RBI = S.statek.RBI;
ra1B = P.sensorParams.raB(:,1);
RpI = P.sensorParams.RpL;
r0G = P.sensorParams.r0G;
ra2B = P.sensorParams.raB(:,2);
sigmab = P.sensorParams.sigmab;

%Obtain the coordinates of the primary antenna in the I frame
rpI = rI + (RBI')*ra1B;
%Change the coordinates to the G frame
RIG = Recef2enu(r0G);
rpG = (RIG')*rpI;
%Obtain noise vector
RpG = (RIG')*RpI*RIG;
wpG = mvnrnd(zeros(3,1), RpG)';
%Obtain "measured" position
rpGtilde = rpG + wpG;

%Obtain the coordinates of the secondary antenna in the I frame
rsI = rI + (RBI')*ra2B;
%Change the coordinates to the G frame
rsG = (RIG')*rsI;
%Obtain the "baseline vector" pointing from the primary antenna to the secondary antenna
rbG = rsG - rpG;
%Obtain RbG error covariance matrix 
rbGnorm = norm(rbG);
rbGunit = rbG/rbGnorm;
epsillon = 10^(-8);
I3 = diag([1 1 1]);
RbG = rbGnorm^2*sigmab^2*(I3 - rbGunit*(rbGunit')) + epsillon*I3;
%Obtain noise vector
wbG = mvnrnd(zeros(3,1), RbG)';
%Obtain "measured" vector
rbGtilde = rbG + wbG;
rbGtildeNorm = norm(rbGtilde);
rbGtilde = rbGnorm*rbGtilde/rbGtildeNorm;

end