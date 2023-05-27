% Top-level script for calling simulateQuadrotorDynamics
clear; clc; close all;
% Quadrotor parameters and constants
quadParamsScript;
constantsScript;
S.quadParams = quadParams;
S.constants = constants;
% Total simulation time, in seconds
Tsim = 4;
% Update interval, in seconds.  This value should be small relative to the
% shortest time constant of your system.
delt = 0.005;
% Time vector, in seconds 
N = floor(Tsim/delt);
S.tVec = [0:N-1]'*delt;
% Matrix of disturbance forces acting on the body, in Newtons, expressed in I
S.distMat = zeros(N-1,3);
% Rotor speeds at each time, in rad/s
S.omegaMat = ones(N-1,4);
S.omegaMat(:,1) =  588.0203*S.omegaMat(:,1);
S.omegaMat(:,2) =  587.8935*S.omegaMat(:,2);
S.omegaMat(:,3) =  587.8935*S.omegaMat(:,3);
S.omegaMat(:,4) =  588.0203*S.omegaMat(:,4);
% Initial position in m
S.state0.r = [0 -0.5 0]';
% Initial attitude expressed as Euler angles, in radians
S.state0.e = [-0.1252 0 0]';
% Initial velocity of body with respect to I, expressed in I, in m/s
S.state0.v = [pi/4 0 0]';
% Initial angular rate of body with respect to I, expressed in B, in rad/s
S.state0.omegaB = [0 -0.1962 1.5585]';
% Oversampling factor
S.oversampFact = 10;
P = simulateQuadrotorDynamics(S);
S2.tVec = P.tVec;
S2.rMat = P.state.rMat;
S2.eMat = P.state.eMat;
S2.plotFrequency = 20;
S2.makeGifFlag = false;
S2.gifFileName = 'testGif.gif';
S2.bounds=1*[-1 1 -1 1 -1 1];
visualizeQuad(S2);
figure(1);clf;
plot(P.tVec,P.state.rMat(:,3)); grid on;
xlabel('Time (sec)');
ylabel('Vertical (m)');
title('Vertical position of CM'); 
figure(2);clf;
plot(P.state.rMat(:,1), P.state.rMat(:,2)); 
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Horizontal position of CM');