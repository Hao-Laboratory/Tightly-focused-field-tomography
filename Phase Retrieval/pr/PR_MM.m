function [amp,phs] = PR_MM(Obj,PSF,Scope,pupilRes,Option)
%PR_INT calculates phase retrieval using integral (scalar version)
%
% *************************************************************************
% Author: Xin Liu
% Emial:liuxin2018@zju.edu.cn
% Jan.07, 2021

if strcmp(Option.Precision,'single')
    Obj.NA = single(Obj.NA);
    Obj.n = single(Obj.n);
    PSF.wavelength = single(PSF.wavelength);
    PSF.amp = single(PSF.amp);
    PSF.phs = single(PSF.phs);
    Scope.xs = single(Scope.xs);
    Scope.ys = single(Scope.ys);
    Scope.zs = single(Scope.zs);
end

%% gpu data preparation
if Option.UseGpu == 1
    Obj.NA = gpuArray(Obj.NA);
    Obj.n = gpuArray(Obj.n);
    PSF.wavelength = gpuArray(PSF.wavelength);
    PSF.amp = gpuArray(PSF.amp);
    PSF.phs = gpuArray(PSF.phs);
    Scope.xs = gpuArray(Scope.xs);
    Scope.ys = gpuArray(Scope.ys);
    Scope.zs = gpuArray(Scope.zs);
end

%% pupil sampling
r0 = 1;  % pupil is defined as an unit circle
[xx,yy] = meshgrid(linspace(-r0,r0,pupilRes));  % cartesian coordinate of pupil plane
[~,rho] = cart2pol(xx,yy);  % polar coordinate of pupil plane

thetaMax = asin(Obj.NA/Obj.n);  % maximum convergent angle of objective
sinTheta = sin(thetaMax).*rho;
sinTheta(rho>1) = 0;
theta = asin(sinTheta);  % convergance angle of each ray(theta)

%% wave vector of time reversal
k0 = 2*pi/PSF.wavelength;  % wavenumber in vacuum
kc = k0*Obj.NA;  % cut-off frequency of objective (in k-space)
kx = linspace(-kc,kc,pupilRes);  % spatial frequency coordinates
ky = kx;
kz = k0*Obj.n*cos(theta);

%% position vector of dipoles
lz = length(Scope.zs);
Scope.zs = reshape(Scope.zs,1,1,lz);

%% electric field distribution of focal space
Ef = PSF.amp.*exp(-1i*PSF.phs);

%% time reversal
Mx = exp(-1i*kx.'*Scope.xs);
My = exp(-1i*Scope.ys.'*ky);
PHI = exp(1i*kz.*Scope.zs);

Ep = pagemtimes(pagemtimes(My.',Ef),Mx.');
% Ep = mdft(Ef,Scope.xs,Scope.ys,kx/(2*pi),ky/(2*pi));

Ep = Ep.*PHI;
Ep = mean(Ep,3);
Ep = conj(Ep);

% apodization for energy conservation
AFF = sqrt(cos(theta)./Obj.n);

% amplitude projection factor from (theta, phi) to (kx, ky) coordinate
APF = 1./kz;

% prefix (the same as that in forward focusing)
f = r0/sin(thetaMax);  % focal length
k = k0*Obj.n;
prefix = -1i*f*exp(1i*k*f)/(2*pi);
dkxdky2dfxdfy = (2*pi)^2;  % convert wave vector to spatial frequency

% pupil function
Ep = Ep./AFF./APF./prefix./dkxdky2dfxdfy;

% energy conservation
dx = (Scope.xs(end) - Scope.xs(1))/(length(Scope.xs)-1);
dy = (Scope.ys(end) - Scope.ys(1))/(length(Scope.ys)-1);
if length(Scope.ys)==1
    dy = 1;
end
if length(Scope.xs)==1
    dx = 1;
end
Ep = Ep.*dx*dy;

Ep(~isfinite(Ep)) = 0;

amp = abs(Ep);
phs = wrapTo2Pi(angle(Ep));

amp(rho>1) = 0;
phs(rho>1) = 0;

if Option.UseGpu == 1
    amp = gather(amp);
    phs = gather(phs);
end

if strcmp(Option.Precision,'single')
    amp = double(amp);
    phs = double(phs);
end

end