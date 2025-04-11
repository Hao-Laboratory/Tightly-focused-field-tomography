%SINGLEOBJECTIVEPSF_NOT_MM calculates the PSF of single objective by
%matrix triple product in Cartesian coordinate system.
%
% INPUT********************************************************************
% Obj.NA: scalar value, numerical aperture of objective lens
% Obj.n: scalar value, refractive index in focal space
% Beam.wavelength: scalar value, Beam.wavelength of light
% Beam.amp: M*M matrix, amplitude distribution on pupil plane
% Beam.phs: M*M matrix, phase distribution on pupil plane
% Beam.abr: M*M matrix, aberration distribution on pupil plane
% Beam.plr: M*M*3 matrix, polarization distribution on pupil plane
% Scope.xs: 1*N array, representing x axis of PSF
% Scope.ys: 1*N array, representing y axis of PSF
% Scope.zs: 1*N array, representing z axis of PSF
% Beam.PupilRes: scalar value, resolution of pupil
% Option.UseGpu: 0 or 1, option using GPU acceleration
% Option.Precision: 0 or 1, precision of numbers
%
% OUTPUT*******************************************************************
% [Ex,Ey,Ez]: electric field of focused field
%
% *************************************************************************
% LIU Xin
% liuxin2018@zju.edu.cn
% Apr.24, 2021

function [Ex,Ey,Ez] = singleobjectivepsf_noT_MM(Obj,Beam,Scope,Option)
%% data initialization
if strcmp(Option.Precision,'single')
    Obj.NA = single(Obj.NA);
    Obj.n = single(Obj.n);
    Beam.wavelength = single(Beam.wavelength);
    Beam.amp = single(Beam.amp);
    Beam.phs = single(Beam.phs);
    Beam.abr = single(Beam.abr);
    Beam.plr = single(Beam.plr);
    Scope.xs = single(Scope.xs);
    Scope.ys = single(Scope.ys);
    Scope.zs = single(Scope.zs);
end

%% gpu data preparation
if Option.UseGpu == 1
    Obj.NA = gpuArray(Obj.NA);
    Obj.n = gpuArray(Obj.n);
    Beam.wavelength = gpuArray(Beam.wavelength);
    Beam.amp = gpuArray(Beam.amp);
    Beam.phs = gpuArray(Beam.phs);
    Beam.abr = gpuArray(Beam.abr);
    Beam.plr = gpuArray(Beam.plr);
    Scope.xs = gpuArray(Scope.xs);
    Scope.ys = gpuArray(Scope.ys);
    Scope.zs = gpuArray(Scope.zs);
end

%% pupil sampling
r0 = 1;  % radius of the pupil is unit
[xx,yy] = meshgrid(linspace(-r0,r0,Beam.PupilRes));  % cartesian coordinate of pupil plane
[phi,rho] = cart2pol(xx,yy);  % polar coordinate of pupil plane

thetaMax = asin(Obj.NA/Obj.n);  % maximum convergent angle of objective
sinTheta = sin(thetaMax).*rho;  % r = f*sin(theta); only for objective obeying sine condition
sinTheta(rho>1) = 0;
theta = asin(sinTheta);  % convergance angle of each ray(theta)

%% interpolation of pupil
if ~(size(Beam.amp) == Beam.PupilRes)
    Beam.amp = pupilInterp(Beam.amp,Beam.PupilRes);
end
if ~(size(Beam.phs) == Beam.PupilRes)
    Beam.phs = pupilInterp(Beam.phs,Beam.PupilRes);
end
if ~(size(Beam.abr) == Beam.PupilRes)
    Beam.abr = pupilInterp(Beam.abr,Beam.PupilRes);
end

px = Beam.plr(:,:,1);
py = Beam.plr(:,:,2);
% pz = Beam.plr(:,:,3);

if ~(size(px) == Beam.PupilRes)
    px = pupilInterp(px,Beam.PupilRes);
end

if ~(size(py) == Beam.PupilRes)
    py = pupilInterp(py,Beam.PupilRes);
end

% if ~(size(pz) == Beam.PupilRes)
%     pz = pupilInterp(pz,Beam.PupilRes);
% end

%% pupil function
E_inc = Beam.amp.*exp(1i.*(Beam.phs+Beam.abr));

% remove parts outside numerical aperture
E_inc(rho>r0) = 0;

fc = Obj.NA/Beam.wavelength;  % cut-off frequency of objective (in k-space)
df = 2*fc/(Beam.PupilRes-1);  % sampling period in k-space

%% polarization in image space (exit pupil)
[Px,Py,Pz] = coortrans(px,py,0,theta,phi,'o2i');

%% amplitude apodization function for energy conservation
% different objective may have different AF!!!
AF = sqrt(cos(theta)./Obj.n);  % only for objective obeying the sine condition

%% frequency domain
k0 = 2*pi/Beam.wavelength;  % wavenumber in vacuum
kc = k0*Obj.NA;  % cut-off frequency of objective (in k-space)

% note that signs of kx and ky are opposite to xp and yp, i.e.,
% kx = -sin(theta).*cos(phi) = -kc*xp/r0; ky = -sin(theta).*sin(phi) = -kc*yp/r0;
kx = linspace(-kc,kc,Beam.PupilRes);  % spatial frequency coordinates
ky = kx;
kz = k0*Obj.n*cos(theta);

%% amplitude transformation
E_inf = E_inc.*AF;

%% the coefficient in front of Debye integral
f = r0/sin(thetaMax);  % focal length (r0 = f*n*sin(thetaMax)?)
k = 2*pi*Obj.n/Beam.wavelength;
prefix = -1i*f*exp(1i*k*f)/(2*pi);
dkxdky2dfxdfy = (2*pi)^2;  % convert wave vector to spatial frequency

% effective scalar field for Fourier transform
E_eff = dkxdky2dfxdfy*prefix*E_inf./kz;

%% plane waves in image space (three polarization components)
Ex0 = E_eff.*Px;
Ey0 = E_eff.*Py;
Ez0 = E_eff.*Pz;

%% initialization
lz = length(Scope.zs);

%% matrix multiplication preparation (dot(k,r))
Mx = exp(-1i*kx.'*Scope.xs);
My = exp(-1i*Scope.ys.'*ky);

%% complexity cost
M = Beam.PupilRes;
N = Beam.PupilRes;
R = length(Scope.ys);
S = length(Scope.xs);

CostMM_Right = M*S*(N+R);
CostMM_Left = N*R*(M+S);

%% batch segmentation
% avoid Option.BatchSize = 0
if Option.BatchSize > lz
    Option.BatchSize = lz;
end

% number of blocks
batchNumber = floor(lz/Option.BatchSize);

% length of the last block
lastBatchSize = mod(lz,Option.BatchSize);

%% creat block matrix
batchVector = ones(1,batchNumber).*Option.BatchSize;
if lastBatchSize ~= 0
    batchVector = [batchVector,lastBatchSize];
end

totalBatches = length(batchVector);  % total cycles

zc = mat2cell(reshape(Scope.zs,1,1,lz),1,1,batchVector);

% memory preallocation
Ex = cell(1,totalBatches);
Ey = Ex;
Ez = Ex;

%% intergral
for ii = 1:totalBatches
    defocusTerm = exp(1i*kz.*zc{ii});
    
    Ewx = Ex0.*defocusTerm;
    Ewy = Ey0.*defocusTerm;
    Ewz = Ez0.*defocusTerm;
    
    if CostMM_Left<=CostMM_Right  % choose the optimal order
        Ex{ii} = pagemtimes(pagemtimes(My,Ewx),Mx);
        Ey{ii} = pagemtimes(pagemtimes(My,Ewy),Mx);
        Ez{ii} = pagemtimes(pagemtimes(My,Ewz),Mx);
    else
        Ex{ii} = pagemtimes(My,pagemtimes(Ewx,Mx));
        Ey{ii} = pagemtimes(My,pagemtimes(Ewy,Mx));
        Ez{ii} = pagemtimes(My,pagemtimes(Ewz,Mx));
    end
end

Ex = cat(3,Ex{:});
Ey = cat(3,Ey{:});
Ez = cat(3,Ez{:});

% energy normalization
Ex = Ex.*df^2;
Ey = Ey.*df^2;
Ez = Ez.*df^2;

if Option.UseGpu == 1
    Ex = gather(Ex);
    Ey = gather(Ey);
    Ez = gather(Ez);
end

if strcmp(Option.Precision,'single')
    Ex = double(Ex);
    Ey = double(Ey);
    Ez = double(Ez);
end
end