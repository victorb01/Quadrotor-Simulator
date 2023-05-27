function Q = simulateQuadrotorControl(R,S,P)
% simulateQuadrotorControl : Simulates closed-loop control of a quadrotor
%                            aircraft.
%
%
% INPUTS
%
% R ---------- Structure with the following elements:
%
%          tVec = Nx1 vector of uniformly-sampled time offsets from the
%                 initial time, in seconds, with tVec(1) = 0.
%
%        rIstar = Nx3 matrix of desired CM positions in the I frame, in
%                 meters.  rIstar(k,:)' is the 3x1 position at time tk =
%                 tVec(k).
%
%        vIstar = Nx3 matrix of desired CM velocities with respect to the I
%                 frame and expressed in the I frame, in meters/sec.
%                 vIstar(k,:)' is the 3x1 velocity at time tk = tVec(k).
%
%        aIstar = Nx3 matrix of desired CM accelerations with respect to the I
%                 frame and expressed in the I frame, in meters/sec^2.
%                 aIstar(k,:)' is the 3x1 acceleration at time tk =
%                 tVec(k).
%
%        xIstar = Nx3 matrix of desired body x-axis direction, expressed as a
%                 unit vector in the I frame. xIstar(k,:)' is the 3x1
%                 direction at time tk = tVec(k).
%  
% S ---------- Structure with the following elements:
%
%  oversampFact = Oversampling factor. Let dtIn = R.tVec(2) - R.tVec(1). Then
%                 the output sample interval will be dtOut =
%                 dtIn/oversampFact. Must satisfy oversampFact >= 1.
%
%        state0 = State of the quad at R.tVec(1) = 0, expressed as a structure
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
%                 disturbance vector acting on the quad from R.tVec(k) to
%                 R.tVec(k+1).
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
% Q ---------- Structure with the following elements:
%
%          tVec = Mx1 vector of output sample time points, in seconds, where
%                 Q.tVec(1) = R.tVec(1), Q.tVec(M) = R.tVec(N), and M =
%                 (N-1)*oversampFact + 1.
%  
%         state = State of the quad at times in tVec, expressed as a
%                 structure with the following elements:
%                   
%                rMat = Mx3 matrix composed such that rMat(k,:)' is the 3x1
%                       position at tVec(k) in the I frame, in meters.
% 
%                eMat = Mx3 matrix composed such that eMat(k,:)' is the 3x1
%                       vector of Euler angles at tVec(k), in radians,
%                       indicating the attitude.
%
%                vMat = Mx3 matrix composed such that vMat(k,:)' is the 3x1
%                       velocity at tVec(k) with respect to the I frame
%                       and expressed in the I frame, in meters per
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

%Find simulation length
tVecSZ = size(R.tVec, 1);

%Time calculations
dtIn = R.tVec(2) - R.tVec(1);
dtOut = dtIn/S.oversampFact;

%Size of output states
finalSZ = (tVecSZ - 1)*S.oversampFact + 1;

%Initialize final output state
Q.state.rMat = zeros(finalSZ, 3);
Q.state.eMat = zeros(finalSZ, 3);
Q.state.vMat = zeros(finalSZ, 3);
Q.state.omegaBMat = zeros(finalSZ, 3);
Q.state.omegaVec = zeros(finalSZ, 4);
Q.tVec = zeros(finalSZ, 1);

%Pack initial state
Q.state.rMat(1,:) = S.state0.r(:);
Q.state.eMat(1,:) = S.state0.e(:);
Q.state.vMat(1,:) = S.state0.v(:);
Q.state.omegaBMat(1,:) = S.state0.omegaB(:);
Q.state.omegaVec(1,:) = zeros(1,4);

Q.tVec(1) = R.tVec(1);

for i = 1:(tVecSZ-1)
    %Define tspan for each iteration
    tspan = (R.tVec(i):dtOut:R.tVec(i+1))';

    %Define step count
    k = (i-1)*S.oversampFact + 1;

    %Unpack distVec for each iteration
    distVec = S.distMat(i,:)';

    %Unpack previous state
    rI = Q.state.rMat(k,:)';
    e = Q.state.eMat(k,:)';
    RBI = euler2dcm(e);
    RBIVec = RBI(:);
    vI = Q.state.vMat(k,:)';
    omegaB = Q.state.omegaBMat(k,:)';
    omegaVec = Q.state.omegaVec(k,:)';

    %Assemble S1 structure
    S1.statek.rI = rI;
    S1.statek.vI = vI;
    S1.statek.RBI = RBI;
    S1.statek.omegaB = omegaB;

    %Assemble R1 structure
    R1.rIstark = R.rIstar(i,:)';
    R1.vIstark = R.vIstar(i,:)';
    R1.aIstark = R.aIstar(i,:)';
    R1.xIstark = R.xIstar(i,:)';

    %Get Fk and zIstark
    [Fk, zIstark] = trajectoryController(R1, S1, P);

    %Get NBk
    R1.zIstark = zIstark;
    NBk = attitudeController(R1, S1, P);

    %Get voltages for each iteration
    eaVec = voltageConverter(Fk, NBk, P);

    %Put together previous state vector
    Xi = [rI; vI; RBIVec; omegaB; omegaVec];

    %Make function handle
    Xdot = @(t, X)quadOdeFunctionHF(t, X, eaVec, distVec, P);

    %Call ode45
    [tVecTemp, stateMatTemp] = ode45(Xdot, tspan, Xi);

    %Save new time variables in P.tVec
    tempSZ = size(tVecTemp);
    Q.tVec(k:(k+tempSZ-1),:) = tVecTemp;

    %Save new states
    rIf = stateMatTemp(:,1:3);
    vIf = stateMatTemp(:,4:6);
    omegaBf = stateMatTemp(:,16:18);
    omegaVecf = stateMatTemp(:,19:22);

    for j = 1:tempSZ
        RBIf = zeros(3,3);
        RBIf(:) = stateMatTemp(j,7:15);
        ef = dcm2euler(RBIf);
        Q.state.eMat(k+j-1,:) = ef';
    end
    Q.state.rMat(k:(k+tempSZ-1),:) = rIf;
    Q.state.vMat(k:(k+tempSZ-1),:) = vIf;
    Q.state.omegaBMat(k:(k+tempSZ-1),:) = omegaBf;
    Q.state.omegaVec(k:(k+tempSZ-1),:) = omegaVecf;

end