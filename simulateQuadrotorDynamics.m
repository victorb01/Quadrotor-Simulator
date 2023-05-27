function P = simulateQuadrotorDynamics(S)
% simulateQuadrotorDynamics : Simulates the dynamics of a quadrotor aircraft.
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
%  
%      omegaMat = (N-1)x4 matrix of rotor speed inputs.  omegaMat(k,j) is the
%                 constant (zero-order-hold) rotor speed setpoint for the jth rotor
%                 over the interval from tVec(k) to tVec(k+1).
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
%
%+-----------------------------------------------------------------------------+
%Create parameters structure A
A.quadParams = S.quadParams;
A.constants = S.constants;

tVecSZ = size(S.tVec);

dtIn = S.tVec(2) - S.tVec(1);
dtOut = dtIn/S.oversampFact;

P.state.rMat(1,:) = S.state0.r(:);
P.state.eMat(1,:) = S.state0.e(:);
P.state.vMat(1,:) = S.state0.v(:);
P.state.omegaBMat(1,:) = S.state0.omegaB(:);

P.tVec(1) = S.tVec(1);

for i = 1:(tVecSZ-1)
    %Define tspan for each iteration
    tspan = [S.tVec(i):dtOut:S.tVec(i+1)]';

    %Define step count
    k = (i-1)*S.oversampFact + 1;

    %Unpack omegaVec and distVec for each iteration
    omegaVec = transpose(S.omegaMat((i),:));
    distVec = transpose(S.distMat((i),:));

    %Unpack previous state
    rI = transpose(P.state.rMat((k),:));
    e = transpose(P.state.eMat((k),:));
    RBI = euler2dcm(e);
    RBIVec = RBI(:);
    vI = transpose(P.state.vMat((k),:));
    omegaB = transpose(P.state.omegaBMat((k),:));

    %Put together previous state vector
    Xi = [rI; vI; RBIVec; omegaB];

    %Make function handle
    Xdot = @(t, X)quadOdeFunction(t, X, omegaVec, distVec, A);

    %Call ode45
    [tVecTemp, stateMatTemp] = ode45(Xdot, tspan, Xi);

    %Save new time variables in P.tVec
    P.tVec(end) = [];
    P.tVec = [P.tVec; tVecTemp];

    %Save new states
    rIf = stateMatTemp(:,1:3);
    vIf = stateMatTemp(:,4:6);
    omegaBf = stateMatTemp(:,16:18);

    tempSZ = size(tVecTemp);
    for j = 1:tempSZ
        RBIf = zeros(3,3);
        RBIf(:) = stateMatTemp(j,7:15);
        ef = dcm2euler(RBIf);
        P.state.eMat(k+j-1,:) = transpose(ef);
    end
    P.state.rMat(end,:) = [];
    P.state.rMat = [P.state.rMat; rIf];
    P.state.vMat(end,:) = [];
    P.state.vMat = [P.state.vMat; vIf];
    P.state.omegaBMat(end,:) = [];
    P.state.omegaBMat = [P.state.omegaBMat; omegaBf];
end