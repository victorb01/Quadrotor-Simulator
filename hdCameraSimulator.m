function [rx] = hdCameraSimulator(rXI,S,P)
% hdCameraSimulator : Simulates feature location measurements from the
%                     quad's high-definition camera. 
%
%
% INPUTS
%
% rXI -------- 3x1 location of a feature point expressed in I in meters.
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
% OUTPUTS
%
% rx --------- 2x1 measured position of the feature point projection on the
%              camera's image plane, in pixels.  If the feature point is not
%              visible to the camera (the ray from the feature to the camera
%              center never intersects the image plane, or the feature is
%              behind the camera), then rx is an empty matrix.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+

%Initialise given variables
rXIhomo = [rXI; 1];
RBI = S.statek.RBI;
rI = S.statek.rI;
RCB = P.sensorParams.RCB;
rocB = P.sensorParams.rocB;
K = P.sensorParams.K;
Rc = P.sensorParams.Rc;
p = P.sensorParams.pixelSize;
planeSize = P.sensorParams.imagePlaneSize;

%Obtain I to C frame rotation
RCI = RCB*RBI;
%Obtain translation from I to C origin
t = -RCB*(RBI*rI + rocB);
%Obtain image homogenous coordinates in C frame
x = K*[RCI, t]*rXIhomo;
%Obtain image coordinates in C frame
xc = [(x(1)/x(3)); (x(2)/x(3))];
%Obtain noise vector
wc = mvnrnd(zeros(2,1), Rc)';
%Obtain measured feature
xcTilde = xc/p + wc;
%Check that measured feature lies on the image plane and is in front of the
%camera
if  (abs(xc(1)) < planeSize(1)/2 && abs(xc(2)) < planeSize(2)/2 && x(3) > 0)
    rx = xcTilde;
else
    rx = [];
end

end