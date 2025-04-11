function phs = vortexPhasePlate(pupilRes,l)
% VORTEXPHASEPLATE: generates a vortex phase plate
%
% -------------------------------------------------------------------------
% MatrixSize: the size of Matrix (must be odd).
% l: the topological charge.
% IsAntiClockWise: the rotation direction. If it is false, the rotation is
% clockwise, otherwise, anticlockwise.
% -------------------------------------------------------------------------
% first coded by Hao,Xiang on Jan.29, 2015
% last updated by Hao,Xiang on Jan.29, 2015
% finally updated by Liu,Xin on Feb.20, 2019
% -------------------------------------------------------------------------

X = linspace(-1,1,pupilRes);
Y = X;
[x,y] = meshgrid(X,Y);
[phi,rho] = cart2pol(x,y);

phs = phi*l;

phs(rho>1) = 0;
phs = mod(phs,2*pi);
end