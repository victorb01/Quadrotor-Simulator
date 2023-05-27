function [P] = simulateQuadrotorDynamicsHF(S)
% simulateQuadrotorDynamicsHF : Simulates the dynamics of a quadrotor
%                               aircraft (high-fidelity version).
%
%
% INPUTS
%
% S ---------- Structure with the following elements:
%
%          tVec = Nx1 vector of uniformly-sampled time offsets from the
%                 initial time, in seconds, with tVec(1) = 0.
%
%  oversampFact = Oversampling factor. Let dtIn = tVec(2) - tVec(1). Then the
%                 output sample interval will be dtOut =
%                 dtIn/oversampFact. Must satisfy oversampFact >= 1.   
%
%         eaMat = (N-1)x4 matrix of motor voltage inputs.  eaMat(k,j) is the
%                 constant (zero-order-hold) voltage for the jth motor over
%                 the interval from tVec(k) to tVec(k+1).
%
%        state0 = State of the quad at tVec(1) = 0, expressed as a structure
%                 with the following elements:
%                   
%                   r = 3x1 position in the world frame, in meters
% 
%                   e = 3x1 vector of Euler angles, in radians, indicating the
%                       attitude
%
%                   v = 3x1 velocity with respect to the world frame and
%                       expressed in the world frame, in meters per second.
%                 
%              omegaB = 3x1 angular rate vector expressed in the body frame,
%                       in radians per second.
%
%       distMat = (N-1)x3 matrix of disturbance forces acting on the quad's
%                 center of mass, expressed in Newtons in the world frame.
%                 distMat(k,:)' is the constant (zero-order-hold) 3x1
%                 disturbance vector acting on the quad from tVec(k) to
%                 tVec(k+1).
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
% P ---------- Structure with the following elements:
%
%          tVec = Mx1 vector of output sample time points, in seconds, where
%                 P.tVec(1) = S.tVec(1), P.tVec(M) = S.tVec(N), and M =
%                 (N-1)*oversampFact + 1.
%                  
%  
%         state = State of the quad at times in tVec, expressed as a structure
%                 with the following elements:
%                   
%                rMat = Mx3 matrix composed such that rMat(k,:)' is the 3x1
%                       position at tVec(k) in the world frame, in meters.
% 
%                eMat = Mx3 matrix composed such that eMat(k,:)' is the 3x1
%                       vector of Euler angles at tVec(k), in radians,
%                       indicating the attitude.
%
%                vMat = Mx3 matrix composed such that vMat(k,:)' is the 3x1
%                       velocity at tVec(k) with respect to the world frame
%                       and expressed in the world frame, in meters per
%                       second.
%                 
%           omegaBMat = Mx3 matrix composed such that omegaBMat(k,:)' is the
%                       3x1 angular rate vector expressed in the body frame in
%                       radians, that applies at tVec(k).
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+  

%Create parameters structure A
A.quadParams = S.quadParams;
A.constants = S.constants;

%Find simulation length
tVecSZ = size(S.tVec, 1);

%Time calculations
dtIn = S.tVec(2) - S.tVec(1);
dtOut = dtIn/S.oversampFact;

%Size of output states
finalSZ = (tVecSZ - 1)*S.oversampFact + 1;

%Initialize final output state
P.state.rMat = zeros(finalSZ, 3);
P.state.eMat = zeros(finalSZ, 3);
P.state.vMat = zeros(finalSZ, 3);
P.state.omegaBMat = zeros(finalSZ, 3);
P.state.omegaVec = zeros(finalSZ, 4);
P.tVec = zeros(finalSZ, 1);

%Pack initial state
P.state.rMat(1,:) = S.state0.r(:);
P.state.eMat(1,:) = S.state0.e(:);
P.state.vMat(1,:) = S.state0.v(:);
P.state.omegaBMat(1,:) = S.state0.omegaB(:);
P.state.omegaVec(1,:) = zeros(1,4);

P.tVec(1) = S.tVec(1);

for i = 1:(tVecSZ-1)
    %Define tspan for each iteration
    tspan = [S.tVec(i):dtOut:S.tVec(i+1)]';

    %Define step count
    k = (i-1)*S.oversampFact + 1;

    %Unpack eaVec and distVec for each iteration
    eaVec = transpose(S.eaMat(i,:));
    distVec = transpose(S.distMat(i,:));

    %Unpack previous state
    rI = transpose(P.state.rMat(k,:));
    e = transpose(P.state.eMat(k,:));
    RBI = euler2dcm(e);
    RBIVec = RBI(:);
    vI = transpose(P.state.vMat(k,:));
    omegaB = transpose(P.state.omegaBMat(k,:));
    omegaVec = transpose(P.state.omegaVec(k,:));

    %Put together previous state vector
    Xi = [rI; vI; RBIVec; omegaB; omegaVec];

    %Make function handle
    Xdot = @(t, X)quadOdeFunctionHF(t, X, eaVec, distVec, A);

    %Call ode45
    [tVecTemp, stateMatTemp] = ode45(Xdot, tspan, Xi);

    %Save new time variables in P.tVec
    tempSZ = size(tVecTemp);
    P.tVec(k:(k+tempSZ-1),:) = tVecTemp;

    %Save new states
    rIf = stateMatTemp(:,1:3);
    vIf = stateMatTemp(:,4:6);
    omegaBf = stateMatTemp(:,16:18);
    omegaVecf = stateMatTemp(:,19:22);

    for j = 1:tempSZ
        RBIf = zeros(3,3);
        RBIf(:) = stateMatTemp(j,7:15);
        ef = dcm2euler(RBIf);
        P.state.eMat(k+j-1,:) = transpose(ef);
    end
    P.state.rMat(k:(k+tempSZ-1),:) = rIf;
    P.state.vMat(k:(k+tempSZ-1),:) = vIf;
    P.state.omegaBMat(k:(k+tempSZ-1),:) = omegaBf;
    P.state.omegaVec(k:(k+tempSZ-1),:) = omegaVecf;
end