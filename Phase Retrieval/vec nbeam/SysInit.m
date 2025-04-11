function [sys, beam] = SysInit(obj,beam,scope,Target,Option)
%SYSINIT implements parameters initialization

if strcmp(Option.Precision,'single')
    obj.NA = single(obj.NA);
    obj.n = single(obj.n);
    beam.wavelength = single(beam.wavelength);
    scope.xs = single(scope.xs);
    scope.ys = single(scope.ys);
    scope.zs = single(scope.zs);
    switch Option.Channel_num
        case 1
            Target.TI = single(Target.TI);
        case 2
            Target.TI_xplr = single(Target.TI_xplr);
            Target.TI_yplr = single(Target.TI_yplr);
    end
end

if Option.UseGpu == 1
    obj.NA = gpuArray(obj.NA);
    obj.n = gpuArray(obj.n);
    beam.wavelength = gpuArray(beam.wavelength);
    scope.xs = gpuArray(scope.xs);
    scope.ys = gpuArray(scope.ys);
    scope.zs = gpuArray(scope.zs);
    
    switch Option.Channel_num
        case 1
            Target.TI = gpuArray(Target.TI);
        case 2
            Target.TI_xplr = gpuArray(Target.TI_xplr);
            Target.TI_yplr = gpuArray(Target.TI_yplr);
    end
end

%% sampling condition of pupil
kc = obj.NA/beam.wavelength;  % cut-off frequency of this system (in k-space)
dk = 2*kc/(beam.PupilRes-1);  % sampling period in k-space

%% generate pupil grid (theta, phi, rho)
r0 = 1;  % radius of pupil
[xx,yy] = meshgrid(linspace(-r0,r0,beam.PupilRes));  % cartesian coordinate of pupil plane
[phi,rho] = cart2pol(xx,yy);  % polar coordinate of pupil plane

thetaMax = asin(obj.NA/obj.n);  % maximum convergent angle of objective
sinTheta = sin(thetaMax)*rho;
sinTheta(rho>r0) = 0;
theta = asin(sinTheta);  % convergance angle of each ray (theta)
k = 2*pi*obj.n./beam.wavelength;
kz = k*cos(theta);

%% amplitude apodization factor for energy conservation
% different objective may have different AF!!!
AAF = sqrt(cos(theta)/obj.n);  % only for objective obeying the sine condition

%% plane waves in image space
% amplitude projection factor from (theta, phi) to (kx, ky) coordinate for
% integral (FFT)
APF = 1./kz;

%% the coefficient in front of Debye integral
f = r0/sin(thetaMax);  % focal length
dkxdky2dfxdfy = (2*pi)^2;  % convert wave vector to spatial frequency
prefix = -1i*f*exp(1i*k*f)/(2*pi)*dkxdky2dfxdfy;

lx = length(scope.xs); ly = length(scope.ys); lz = length(scope.zs);

% memory preallocation for retrieved PSF
PSFr = zeros(ly,lx,lz);
if strcmp(Option.Precision,'single')
    PSFr = single(PSFr);
end

if Option.UseGpu == 1
    PSFr = gpuArray(PSFr);
end

%% phase retrieval parameters

% sequence of kx and ky
kxy = linspace(-kc,kc,beam.PupilRes);

%% Mx and My for MTP
sys.Mx = exp(-1i*2*pi*kxy.'*scope.xs);
sys.My = exp(-1i*2*pi*scope.ys.'*kxy);

%% OUTPUT
% common
sys.obj = obj;
sys.scope = scope;
sys.wavelength = beam.wavelength;
sys.kc = kc;
sys.df = dk;
sys.dx = (scope.xs(end) - scope.xs(1))/(lx-1);
sys.dy = (scope.ys(end) - scope.ys(1))/(ly-1);

% for forward focusing
sys.theta = theta;
sys.phi = phi;
sys.prefix = prefix.*AAF.*APF;
sys.lz = lz;
sys.lx = lx;
sys.ly = ly;

% for backward retrieval
switch Option.Channel_num
    case 1
        sys.Target.TI = Target.TI;
    case 2
        sys.Target.TI_xplr = Target.TI_xplr;
        sys.Target.TI_yplr = Target.TI_yplr;
end
sys.pupilRes = beam.PupilRes;
sys.rho = rho;
sys.Ir = PSFr;

% other
sys.theta = theta;
end