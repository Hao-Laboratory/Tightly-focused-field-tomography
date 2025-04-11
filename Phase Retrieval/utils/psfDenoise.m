function [psfDenoised, snr] = psfDenoise(psfStack,varargin)
%PSFDENOISE implements PSF denoising
%   
% Required inputs**********************************************************
% psfStack: input PSF; Nx*Ny*Nz matrix
%
% Optional inputs**********************************************************
% sigRad: radius of signal region
% kernelSizeRatio: filter size
% kernelWidthRatio: filter width
% NnoiseIter: iteration number
%
% Reference: 
% Kromann, Emil B., et al. "Quantitative pupil analysis in stimulated
% emission depletion microscopy using phase retrieval." Optics letters
% 37.11 (2012): 1805-1807.
%
% LIU Xin
% liuxin2018@zju.edu.cn
% Jan. 05, 2024

%% defult input
psfSize = size(psfStack);
if psfSize(1) ~= psfSize(2)
    warning('Input PSF must be square shape.');
end
psfRes = psfSize(1);

try
    psfNum = psfSize(3);
catch
    psfNum = 1;
end

% parse arguments
sigRad = 0.8;  % signal radius
NnoiseIter = 100;

NA = 1.35;
wavelength = 775;
pixelsize = 30;
% th = 1e-2;
th = 0.2;

for ii = 1:2:length(varargin)-1
    switch lower(varargin{ii})
        case 'sigrad'
            sigRad = varargin{ii+1};
        case 'nnoiseiter'
            NnoiseIter = varargin{ii+1};
        case 'na'
            NA = varargin{ii+1};
        case 'wavelength'
            wavelength = varargin{ii+1};
        case 'pixelsize'
            pixelsize = varargin{ii+1};
        case 'threshold'
            th = varargin{ii+1};
        otherwise
            error(['Invalid parameter: ',varargin{ii}]);
    end
end

%% background subtraction
% estimate and subtract noisefloor based on pixels outside the given radius
xys = linspace(-1,1,psfRes);
[xx,yy] = meshgrid(xys);
% mask = (sqrt(xx.^2+yy.^2) >= sigRad);  % identify pixels outside circle
mask = (abs(xx) >= sigRad | abs(yy) >= sigRad);  % identify pixels outside square

% [psfStack, ~] = otf_filter(psfStack,NA,wavelength,pixelsize);
psfNoiseSub = zeros(psfRes,psfRes,psfNum);
bkg_level = zeros(1,1,psfNum);
noise_level = zeros(1,1,psfNum);
signal_level = zeros(1,1,psfNum);
for ii = 1:psfNum  % subtract background noise in all images
    psfTemp = psfStack(:,:,ii);  % extract current image
    psfBkgd = psfTemp(mask);  % identify 'noisy background'
    
    bkg_level(:,:,ii) = mean(psfBkgd);  % estimate background
    noise_level(:,:,ii) = std(psfBkgd);  % estimate noise
    
    % subtract mean of noisy signal, may produce negative values
    psf_no_bg = psfStack(:,:,ii) - bkg_level(:,:,ii);
    psfNoiseSub(:,:,ii) = psf_no_bg;

    psf_no_bg_norm = psf_no_bg./max(psf_no_bg,[],'all');
    signal_region = psf_no_bg_norm>0.5;
    signal_level(:,:,ii) = mean(psf_no_bg(signal_region));
end
snr = signal_level./noise_level;

%% remove margin of signal region, e.g., scanning mirror offset error
% multiply by cleanup function (big flat-topped gaussian-like function)
kernel_fwhm = 0.25*psfRes*(1-sigRad);  % [pixels]
kernel_std = kernel_fwhm/2.355;
G = double(~mask);  % cleanup mask (signal region)

% smoothed cleanup
cleanupFilter = imgaussfilt(G,kernel_std);

% figure('Name','Cleanup Filter');
% subplot(1,2,1);imagesc(cleanupFilter);axis image;colormap(gray);colorbar;
% subplot(1,2,2);plot(cleanupFilter(round(psfRes/2),:));axis tight;

% multiply cleanup filter onto image
psfClean = psfNoiseSub.*cleanupFilter;  % may be unnecessary

%% apply iterative noise suppression (smooth negative pixels)
% define energy neutral filter kernel (pixel sum = 0)
F_eNeu = [0.1 0.15 0.1; 0.15 -1 0.15; 0.1 0.15 0.1];

% figure('Name','Neutral Filter');
% imagesc(F_eNeu);axis image off;colormap(gray);colorbar;

psf_neg_rem = zeros(psfRes,psfRes,psfNum);
for ii = 1:psfNum
    psfTemp = psfClean(:,:,ii);  % extract current image
    for nnoiseIter = 1:NnoiseIter  % run all iterations in noise suppression scheme
        negMask = psfTemp < 0;  % identify pixels with negative values
        negData = psfTemp.*double(negMask);  % isolate negative pixels (positive are set to zero)
        filData = imfilter(negData,F_eNeu,'circular','same','conv');  % apply filter to negative pixels
        psfTemp = psfTemp + filData;  % add smoothed negative pixels (pixel sum  = 0) to original        
    end  % end iterative noise scheme for current image
    psf_neg_rem(:,:,ii) = psfTemp;  % update current image in image stack
end

%% remove high-frequency artifacts
diff_lim_fwhm = wavelength/(2*NA);  % FWHM of PSF
diff_lim_px = diff_lim_fwhm/pixelsize + 1;
kernel_std = diff_lim_px/2.355;
psfDenoised = zeros(psfRes,psfRes,psfNum);
for ii = 1:psfNum
    psf_tmp = psf_neg_rem(:,:,ii);

    % generate a continuous and closed region without holes
    mask_temp = double(psf_tmp > noise_level(:,:,ii));
    mask_filt = imgaussfilt(mask_temp,kernel_std);
    mask_final = double(mask_filt>th);
    psfDenoised(:,:,ii) = psf_tmp.*mask_final;

end

psfDenoised(psfDenoised<0) = 0;
end