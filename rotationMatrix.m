function R = rotationMatrix(a, phi)
% rotationMatrix : Generates the rotation matrix R corresponding to a rotation
% through an angle phi about the axis defined by the unit
% vector aHat. This is a straightforward implementation of
% Euler’s formula for a rotation matrix.
%
% INPUTS
%
% aHat ------- 3-by-1 unit vector constituting the axis of rotation.
%
% phi -------- Angle of rotation, in radians.
%
%
% OUTPUTS
%
% R ---------- 3-by-3 rotation matrix
%
%+------------------------------------------------------------------------------+
    %Creating identity matrix
    I3 = [1 0 0; 0 1 0; 0 0 1];
    %Creating a matrix by multiplying axis of rotation vector by its
    %transpose
    A = a*transpose(a);
    %Final rotation matrix using Euler's formula
    R = cos(phi)*I3 + (1-cos(phi))*A - sin(phi)*crossProductEquivalent(a);
end