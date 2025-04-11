% Xin Liu, et al. "In situ fully vectorial tomography and pupil function retrieval of
% tightly focused fields" Nature Communications xx.xx (2025): xxx.
%
% If you use this code and find our work valuable, please cite our paper.

clear;clc;
close all;

%% system parameters
global_fig_flag = 1;

Option.UseGpu = 1;
Option.Precision = 'single';
Option.BatchSize = 1;
Option.PrBatchSize = 1;
if global_fig_flag == 1
    Option.PltFig = 1;
else
    Option.PltFig = 0;
end

Option.Channel_num = 2;
% Threshold.MaxIter = 1e4; Option.step_size = 1; Option.lr_w = 1e-3;  % optimizable weight
Threshold.MaxIter = 3e3; Option.step_size = 1; Option.lr_w = 0;  % fixed weight

%% colormap
cmap_A = parula;
cmap_P = hsv;
Option.PlrCmap = turbo;

cmap_PSF = hot;  % PSF amplitude
cmap_PSFP = hsv;  % PSF phase

Option.PlrNum = 15;
Option.PlrStyle = 'patch';  % preferred
Option.PlrBgcolor = 'k';  % preferred
Option.PlrLw = 0.2;
Option.PlrNorm = 1;

%% computation parameters
Beam1.PupilRes = 128;
Beam2.PupilRes = 128;
Beam3.PupilRes = 256;

psfRes = 64;
pixelSize = 120;

n_beam = 1;

roi_ratio = 0.6;
rand_Level = 5;
rand_N = 10*n_beam;  % number of random phase modulation

photon_num_peak = 300;
photon_num_bg = 10;
photon_num_add = 0;

%% PR mode
Option.PrMode = ['amp','phs','plr'];

%% objective
Obj.NA = 1.35; Obj.n = 1.518;  % oil sample

% direct detection
Option.w_x = 0.39; Option.w_y = 0.39; Option.w_z = 0.22;

% two-channel (theoretical data)
w = [0.38 0.01 0.11; 0.01 0.38 0.11];

Option.w_x_xplr = w(1,1);
Option.w_y_xplr = w(1,2);
Option.w_z_xplr = w(1,3);

Option.w_x_yplr = w(2,1);
Option.w_y_yplr = w(2,2);
Option.w_z_yplr = w(2,3);

%% wavelength
Beam1.wavelength = 775;
Beam2.wavelength = 590;
Beam3.wavelength = 647;

%% amplitude
peak_int_ratio = [0.4 1 0.6];

Beam1.amp = LaguerreGaussAmp(Beam1.PupilRes,1,1);
Beam2.amp = LaguerreGaussAmp(Beam2.PupilRes,1,1);
Beam3.amp = LaguerreGaussAmp(Beam3.PupilRes,1,1);

Beam1.amp = circularpupil(Beam1.amp);
Beam2.amp = circularpupil(Beam2.amp);
Beam3.amp = circularpupil(Beam3.amp);

%% phase
Beam1.phs = vortexPhasePlate(Beam1.PupilRes,0);
Beam2.phs = vortexPhasePlate(Beam2.PupilRes,0);
Beam3.phs = zeros(Beam3.PupilRes);

%% RMS aberration
Beam1.abr = zeros(Beam1.PupilRes);
Beam2.abr = zeros(Beam2.PupilRes);
Beam3.abr = zeros(Beam3.PupilRes);

%% poalrization
% Beam1.plr = randomPolarization(Beam1.PupilRes,4);  cylinderPlr = 0;
% Beam1.plr = CircularPolarization(Beam1.PupilRes,'l');  cylinderPlr = 0;
Beam1.plr = cylindVecPolarization(Beam1.PupilRes,0,1);  cylinderPlr = 1;
Beam2.plr = cylindVecPolarization(Beam2.PupilRes,pi/2,1);
Beam3.plr = CircularPolarization(Beam3.PupilRes,'l');

%% combine beams
Beams = {Beam1, Beam2, Beam3};
Beam = Beams(1:n_beam);

%% define scope
lz = 5;  % number of z slice
large_roi = (psfRes-1)*pixelSize;
Scope.xs = linspace(-(psfRes-1)/2,(psfRes-1)/2,psfRes)*pixelSize;
Scope.ys = Scope.xs;
Scope.zs = linspace(-1e3,1e3,lz);

dz = linspace(Scope.zs(1),Scope.zs(end),rand_N);  % step for defocus diversity

if global_fig_flag == 1
    % current scope
    fprintf('The current scope is %.2f um. \n', (psfRes-1)*pixelSize*1e-3);
end

%% peak intensity control
% to prevent low signals from being overwhelmed by noise
% use random phase diversity to precalculate intensity of PSFs
Scope_tmp = Scope;
Scope_tmp.zs = 0;

PSF_fp0 = cell(1,n_beam);
for jj = 1:rand_N
    for ii = 1:n_beam
        roi_size_x_tmp = large_roi/2.*roi_ratio;
        roi_size_y_tmp = large_roi/2.*roi_ratio;
        phs_div_tmp = rand_phase_gen_inv_grad(Obj,Beam{ii},roi_size_x_tmp,roi_size_y_tmp,rand_Level);
        Beam_rand = Beam{ii};
        Beam_rand.phs = Beam_rand.phs + phs_div_tmp;
        [Ex,Ey,Ez] = singleobjectivepsf_noT_MM(Obj,Beam_rand,Scope_tmp,Option);
        PSF_fp0{ii}(:,:,jj) = abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;
    end
end

% correct pupil amplitude to control peak intensity ratio among beams
peak_int_ratio_real = zeros(1,n_beam);
for ii = 1:n_beam
    peak_int_ratio_real(ii) = mean(max(PSF_fp0{ii},[],[1,2]));
end

for ii = 1:n_beam
    Beam{ii}.amp = Beam{ii}.amp/sqrt(peak_int_ratio_real(ii))*sqrt(peak_int_ratio(ii));
end

ampTotal = [];
for ii = 1:n_beam
    ampTotal = cat(1,ampTotal,Beam{ii}.amp(:));
end
for ii = 1:n_beam
    Beam{ii}.amp = Beam{ii}.amp./max(ampTotal(:));
end

%% show preset pupil function
if global_fig_flag == 1
    fig_pupil = figure('Name','pupil function');
    t_pupil = tiledlayout(fig_pupil,2*n_beam,3,'TileSpacing','compact','Padding','compact');
    for ii = 1:n_beam
        nexttile(t_pupil,1+(ii-1)*6);
        pupilshow(Beam{ii}.amp);
        colormap(gca,cmap_A);
        colorbar;
        clim([0 1]);
        title('amp GT');

        nexttile(t_pupil,2+(ii-1)*6);
        pupilshow(wrapTo2Pi(Beam{ii}.phs+Beam{ii}.abr));
        colormap(gca,cmap_P);
        clim([0 2*pi]);
        colorbar;
        title('phs GT');

        % nexttile(t_pupil,3+(ii-1)*6);
        % pupilshowplr([],Beam{ii}.plr.*Beam{ii}.amp,'num',Option.PlrNum,'style',Option.PlrStyle,...
        %     'bgcolor',Option.PlrBgcolor,'cmap',Option.PlrCmap,'lw',Option.PlrLw,'normflag',Option.PlrNorm);
        % colorbar;
        % title('plr GT');
    end
end

%% recalculate PSF
I = zeros(psfRes,psfRes,lz);
PSFs = cell(1,n_beam);
PSF_fp = zeros(psfRes,psfRes,n_beam);
for ii = 1:n_beam
    [Ex,Ey,Ez] = singleobjectivepsf_noT_MM(Obj,Beam{ii},Scope,Option);
    PSFs{ii} = abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;
    PSF_fp(:,:,ii) = PSFs{ii}(:,:,round(lz/2));
    I = I + PSFs{ii};
end
PSF = I./max(I(:));

if global_fig_flag == 1
    fig_psfComp = figure('Name','PSF comparison');
    t_psfComp = tiledlayout(fig_psfComp,2,lz,'TileSpacing','compact','Padding','compact');
    for ii = 1:lz
        nexttile(t_psfComp,ii);
        imagesc(Scope.xs,Scope.ys,PSF(:,:,ii));
        axis image xy off;
        colormap(gca,cmap_PSF);
        clim([0 1]);
        title(['z = ',num2str(Scope.zs(ii))]);
        if ii == 1
            axis on;
            xticks([]);
            yticks([]);
            ylabel('PSF GT');
        end
    end
end

%% random phase modulation
ScopeRand = Scope;
ScopeRand.zs = 0;
SLM_phs = cell(1,n_beam);
zeros(Beam1.PupilRes,Beam1.PupilRes,rand_N);

Option.roi_size_x = cell(1,n_beam); Option.roi_size_y = cell(1,n_beam);
tip_tilt = cell(1,n_beam);
Option.x_size = cell(1,n_beam);
Option.y_size = cell(1,n_beam);
for jj = 1:n_beam
    psf_tmp = PSFs{jj}(:,:,round(lz/2));
    [x_cent,y_cent,x_size,y_size] = ...
        psf_cent_size(psf_tmp,0.01,pixelSize,'intuitive');
    roi_size_x = roi_ratio*(large_roi/2);
    roi_size_y = roi_ratio*(large_roi/2);

    Option.roi_size_x{jj} = max([roi_size_x - x_size/2, 0.3*large_roi/2]);
    Option.roi_size_y{jj} = max([roi_size_y - y_size/2, 0.3*large_roi/2]);

    tip_tilt{jj} = focusShiftPhasePlate(Beam{jj}.PupilRes,Beam{jj}.wavelength,Obj.NA,Obj.n,-x_cent,-y_cent,0);
    Option.x_size{jj} = x_size;
    Option.y_size{jj} = y_size;
end

I_f = zeros(psfRes,psfRes,rand_N);
It = zeros(psfRes,psfRes,rand_N);
It_x = zeros(psfRes,psfRes,rand_N);
It_y = zeros(psfRes,psfRes,rand_N);
a_tmp = cell(n_beam,rand_N);
for ii = 1:rand_N
    for jj = 1:n_beam
        rand_phs = rand_phase_gen_inv_grad(Obj,Beam{jj},Option.roi_size_x{jj},Option.roi_size_y{jj},rand_Level);
        rand_phs_tmp = rand_phs + tip_tilt{jj};
        [a_tmp{jj,ii},phs_div_tmp2] = zernikeDecomposition(rand_phs_tmp,50,'norm');
        SLM_phs{jj}(:,:,ii) = phs_div_tmp2;

        Beam_rand = Beam{jj};
        Beam_rand.phs = Beam_rand.phs + SLM_phs{jj}(:,:,ii);
        [Ex,Ey,Ez] = singleobjectivepsf_noT_MM(Obj,Beam_rand,ScopeRand,Option);

        I_f(:,:,ii) = I_f(:,:,ii) + abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;
        It(:,:,ii) = It(:,:,ii) + Option.w_x*abs(Ex).^2 + Option.w_y*abs(Ey).^2 + Option.w_z*abs(Ez).^2;

        It_x(:,:,ii) = It_x(:,:,ii) + Option.w_x_xplr*abs(Ex).^2 + Option.w_y_xplr*abs(Ey).^2 + Option.w_z_xplr*abs(Ez).^2;
        It_y(:,:,ii) = It_y(:,:,ii) + Option.w_x_yplr*abs(Ex).^2 + Option.w_y_yplr*abs(Ey).^2 + Option.w_z_yplr*abs(Ez).^2;
    end
end

It_norm = It./max(It,[],'all');

Itxy = cat(3,It_x,It_y);
Itx_norm = It_x./max(Itxy,[],'all');
Ity_norm = It_y./max(Itxy,[],'all');
Itxy_norm = cat(3,Itx_norm,Ity_norm);

% add noise and denoise
avg_peak_photons_1 = mean(max(It_norm,[],[1,2]));
prefix_1 = photon_num_peak/avg_peak_photons_1;

avg_peak_photons_2 = mean(max(Itxy_norm,[],[1,2]));
prefix_2 = photon_num_peak/avg_peak_photons_2;

% add noise and denoise
It_noise = add_noise(It_norm*prefix_1+photon_num_bg,photon_num_add);
It_noise = It_noise./max(It_noise,[],'all');
It_Denoise = psfDenoise(It_noise,'NA',Obj.NA,'wavelength',Beam1.wavelength,'pixelSize',pixelSize);
It_Denoise = It_Denoise./max(It_Denoise,[],'all');

Itx_noise = add_noise(Itx_norm*prefix_2+photon_num_bg,photon_num_add);
Ity_noise = add_noise(Ity_norm*prefix_2+photon_num_bg,photon_num_add);
Itxy_noise = cat(3,Itx_noise,Ity_noise);
Itx_noise = Itx_noise./max(Itxy_noise,[],'all');
Ity_noise = Ity_noise./max(Itxy_noise,[],'all');
Itx_Denoise = psfDenoise(Itx_noise,'NA',Obj.NA,'wavelength',Beam1.wavelength,'pixelSize',pixelSize);
Ity_Denoise = psfDenoise(Ity_noise,'NA',Obj.NA,'wavelength',Beam1.wavelength,'pixelSize',pixelSize);

% normalization
Itxy_Denoise = cat(3,Itx_Denoise,Ity_Denoise);
Itx_Denoise = Itx_Denoise./max(Itxy_Denoise,[],'all');
Ity_Denoise = Ity_Denoise./max(Itxy_Denoise,[],'all');

if global_fig_flag == 1
    switch Option.Channel_num
        case 1
            I_norm_show = It_norm;
            It_noise_show = It_noise;
            It_Denoise_show = It_Denoise;
        case 2
            I_norm_show = Itx_norm;
            It_noise_show = Itx_noise;
            It_Denoise_show = Itx_Denoise;
    end

    % [~,idx_max] = max(max(I_norm_show,[],[1,2]));
    idx_max = 1;
    figure('Name','sample ideal PSF');
    imagesc(I_norm_show(:,:,idx_max));axis image xy off;colormap(gca,cmap_PSF);title('sample ideal PSF');
    % colormap(gca,parula);
    figure('Name','sample noisy PSF');
    imagesc(It_noise_show(:,:,idx_max));axis image xy off;colormap(gca,cmap_PSF);title('sample noisy PSF');
    % colormap(gca,parula);
    figure('Name','sample denoised PSF');
    imagesc(It_Denoise_show(:,:,idx_max));axis image xy off;colormap(gca,cmap_PSF);title('sample denoised PSF');
    % colormap(gca,parula);

    I_max = max(cat(3,Itx_norm(:,:,idx_max),Ity_norm(:,:,idx_max)),[],'all');
    figure('Name','sample ideal PSF in paper');
    subplot(1,2,1);
    imagesc(Itx_norm(:,:,idx_max));axis image xy off;colormap(gca,cmap_PSF);title('sample ideal PSF');clim([0 I_max]);
    subplot(1,2,2);
    imagesc(Ity_norm(:,:,idx_max));axis image xy off;colormap(gca,cmap_PSF);title('sample ideal PSF');clim([0 I_max]);
end

% energy normalize again (by sum to eliminate power fluctuation of laser)
Target.TI = It_Denoise;
Target.TI_xplr = Itx_Denoise;
Target.TI_yplr = Ity_Denoise;

if global_fig_flag == 1
    N_show = min(11,rand_N);
    fig_randPhsMod = figure('Name','random phase modulation');
    t_randPhsMod = tiledlayout(fig_randPhsMod,4+n_beam,N_show,'TileSpacing','compact','Padding','compact');
    for jj = 1:n_beam
        for ii = 1:N_show
            nexttile(t_randPhsMod,ii + (jj-1)*N_show);
            % pupilshow(SLM_phs{jj}(:,:,ii)); clim([0 max(SLM_phs{jj}(:))]);
            pupilshow(wrapTo2Pi(SLM_phs{jj}(:,:,ii))); clim([0 2*pi]);
            colormap(gca,cmap_P);
            if ii == 1
                axis on;
                xticks([]);
                yticks([]);
                ylabel('random phase');
            end
            title(['P ',num2str(ii)]);
        end
    end

    for ii = 1:N_show
        nexttile(t_randPhsMod,ii+n_beam*N_show);
        imagesc(ScopeRand.xs,ScopeRand.ys,I_norm_show(:,:,ii));
        axis image xy off;
        colormap(gca,cmap_PSF);
        clim([0 1]);
        if ii == 1
            axis on;
            xticks([]);
            yticks([]);
            ylabel('PSF GT');
        end
    end

    for ii = 1:N_show
        nexttile(t_randPhsMod,ii+(n_beam+1)*N_show);
        imagesc(ScopeRand.xs,ScopeRand.ys,It_noise_show(:,:,ii));
        axis image xy off;
        colormap(gca,cmap_PSF);
        clim([0 1]);
        if ii == 1
            axis on;
            xticks([]);
            yticks([]);
            ylabel('PSF noise');
        end
    end

    for ii = 1:N_show
        nexttile(t_randPhsMod,ii+(n_beam+2)*N_show);
        imagesc(ScopeRand.xs,ScopeRand.ys,It_Denoise_show(:,:,ii));
        axis image xy off;
        colormap(gca,cmap_PSF);
        clim([0 1]);
        if ii == 1
            axis on;
            xticks([]);
            yticks([]);
            ylabel('PSF denoised');
        end
    end
end

%% vectorial phase retrieval
ScopeRand.zs = zeros(1,rand_N);
[amp,phs,plr,PSFr, metrics_val, sys, Option] = vec_phasediversity_nbeam(Beam,Obj,ScopeRand,Target,SLM_phs,...
    Threshold,Option);

%% show results
plr_out0 = cell(1,n_beam);
plr_out1 = cell(1,n_beam);

for ii = 1:n_beam
    if contains(Option.PrMode,'phs')
        % phase and polarization visualization correction
        Ep_tmp = plr{ii}.*exp(1i*phs{ii});
        Ep0 = Beam{ii}.plr.*exp(1i*Beam{ii}.phs);
        Ep_crt = complex_correction(Ep_tmp,Ep0);
        [~,phs{ii},plr{ii}] = pupilfuntion_decode(Ep_crt(:,:,1),Ep_crt(:,:,2));
        plr{ii}(:,:,3) = zeros(Beam{ii}.PupilRes);

        [phs{ii},plr{ii}] = phs_plr_correction(Beam{ii}.plr,plr{ii},phs{ii});
        [~,idx_amp_max] = max(amp{ii},[],'all','linear');
        delta_phs_pupil_2d = wrapTo2Pi(phs{ii} - Beam{ii}.phs - Beam{ii}.abr);
        delta_phs_pupil = delta_phs_pupil_2d(idx_amp_max);
        phs{ii} = phs{ii}-delta_phs_pupil;

        plr_out0{ii} = Beam{ii}.plr.*exp(1i*(Beam{ii}.phs + Beam{ii}.abr));
        plr_out1{ii} = plr{ii}.*exp(1i*phs{ii});

        if cylinderPlr == 1  % better visualization for cylindrical polarization
            plr_out0{ii} = plrBasisTransform(plr_out0{ii},'xy2cylinder');
            plr_out1{ii} = plrBasisTransform(plr_out1{ii},'xy2cylinder');
        end
    end
end

if global_fig_flag == 1
    figure('Name','polarization basis comparison');
    tiledlayout(2*n_beam,4,"TileSpacing","compact","Padding","compact");
    for ii = 1:n_beam
        nexttile(1+(ii-1)*8);
        pupilshow(abs(plr_out0{ii}(:,:,1)));
        clim([0 1]);colormap(gca,cmap_A);
        colorbar;title('|E1_{gt}|');
        nexttile(2+(ii-1)*8);
        pupilshow(angle2pi(plr_out0{ii}(:,:,1)));
        clim([0 2*pi]);colorbar;colormap(gca,cmap_P);title('arg(E1_{gt})');
        nexttile(3+(ii-1)*8);
        pupilshow(abs(plr_out0{ii}(:,:,2)));
        clim([0 1]);colormap(gca,cmap_A);
        colorbar;title('|E2_{gt}|');
        nexttile(4+(ii-1)*8);
        pupilshow(angle2pi(plr_out0{ii}(:,:,2)));
        clim([0 2*pi]);colorbar;colormap(gca,cmap_P);title('arg(E2_{gt})');

        nexttile(5+(ii-1)*8);
        pupilshow(abs(plr_out1{ii}(:,:,1)));
        clim([0 1]);
        colorbar;colormap(gca,cmap_A);title('|E1_{ret}|');
        nexttile(6+(ii-1)*8);
        pupilshow(angle2pi(plr_out1{ii}(:,:,1)));
        clim([0 2*pi]);colorbar;colormap(gca,cmap_P);title('arg(E1_{ret})');
        nexttile(7+(ii-1)*8);
        pupilshow(abs(plr_out1{ii}(:,:,2)));
        clim([0 1]);colormap(gca,cmap_A);title('|E2_{ret}|');
        colorbar;
        nexttile(8+(ii-1)*8);
        pupilshow(angle2pi(plr_out1{ii}(:,:,2)));
        clim([0 2*pi]);colorbar;colormap(gca,cmap_P);title('arg(E2_{ret})');
    end
end

%% show retrieved pupil function
if global_fig_flag == 1
    for ii = 1:n_beam
        nexttile(t_pupil,4+(ii-1)*6);
        pupilshow(amp{ii});
        clim([0 1]);
        colormap(gca,cmap_A);
        colorbar;
        title('retrieved amp');

        nexttile(t_pupil,5+(ii-1)*6);
        pupilshow(wrapTo2Pi(phs{ii}));
        colormap(gca,cmap_P);
        clim([0 2*pi]);
        colorbar;
        title('retrieved phs');

        % nexttile(t_pupil,6+(ii-1)*6);
        % pupilshowplr([],plr{ii}.*amp{ii},'num',Option.PlrNum,'style',Option.PlrStyle,...
        %     'bgcolor',Option.PlrBgcolor,'cmap',Option.PlrCmap,'lw',Option.PlrLw,'normflag',Option.PlrNorm);
        % colorbar;
        % title('retrieved plr');
    end
end

%% show reconstructed random modulated PSF
if global_fig_flag == 1
    for ii = 1:N_show
        nexttile(t_randPhsMod,ii+4*N_show);
        imagesc(ScopeRand.xs,ScopeRand.ys,PSFr(:,:,ii)./max(PSFr(:)));
        axis image xy off;
        clim([0 1]);
        colormap(gca,cmap_PSF);
        if ii == 1
            axis on;
            xticks([]);
            yticks([]);
            ylabel('retrieved PSF');
        end
    end

    figure('Name','sample retrieved PSF with modulation');
    imagesc(PSFr(:,:,idx_max));axis image xy off;colormap(gca,cmap_PSF);title('retrieved PSF');
    colormap(gca,cmap_PSF);
end

%% retrieved PSF without random modulation
I_ret = zeros(psfRes,psfRes,lz);
Beam_ret = cell(1,n_beam);
for ii = 1:n_beam
    Beam_ret{ii} = Beam{ii};
    Beam_ret{ii}.phs = phs{ii};
    Beam_ret{ii}.abr = zeros(Beam{ii}.PupilRes);
    Beam_ret{ii}.amp = amp{ii};
    Beam_ret{ii}.plr = plr{ii};

    [Ex,Ey,Ez] = singleobjectivepsf_noT_MM(Obj,Beam_ret{ii},Scope,Option);
    I_ret = I_ret + abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;
end
I_ret = I_ret./max(I_ret(:));

if global_fig_flag == 1
    for ii = 1:lz
        nexttile(t_psfComp,ii+lz);
        imagesc(Scope.xs,Scope.ys,I_ret(:,:,ii));
        axis image xy off;
        clim([0 1]);
        colormap(gca,cmap_PSF);
        if ii == 1
            axis on;
            xticks([]);
            yticks([]);
            ylabel('PSF retrieved');
        end
    end
end

%% show electric components with high resolution in ROI
% define fine grid
small_roi = large_roi/2;

% psfResHi = 256;
psfResHi = 128;
pixelSizeHi = small_roi/(psfResHi-1);

ScopeHi.xs = linspace(-(psfResHi-1)/2,(psfResHi-1)/2,psfResHi)*pixelSizeHi;
ScopeHi.ys = ScopeHi.xs;
ScopeHi.zs = linspace(-1e3,1e3,lz);
% ScopeHi.zs = ScopeHi.xs; lz = psfResHi;

% gt PSF
Ex_gt0 = cell(1,n_beam);
Ey_gt0 = cell(1,n_beam);
Ez_gt0 = cell(1,n_beam);
E_max0 = cell(1,n_beam);
I_gt_all = cell(1,n_beam);
I_gt0 = zeros(psfResHi,psfResHi,lz);
for ii = 1:n_beam
    [Ex_gt0{ii},Ey_gt0{ii},Ez_gt0{ii}] = singleobjectivepsf_noT_MM(Obj,Beam{ii},ScopeHi,Option);
    E_max0{ii} = max(abs(cat(3,Ex_gt0{ii},Ey_gt0{ii},Ez_gt0{ii})),[],"all");
    I_gt_all{ii} = abs(Ex_gt0{ii}).^2 + abs(Ey_gt0{ii}).^2 + abs(Ez_gt0{ii}).^2;
    I_gt0 = I_gt0 + I_gt_all{ii};
    I_gt_all{ii} = I_gt_all{ii}./max(I_gt_all{ii}(:));
end

% normalization
Ex_gt = cell(1,n_beam);
Ey_gt = cell(1,n_beam);
Ez_gt = cell(1,n_beam);
for ii = 1:n_beam
    Ex_gt{ii} = Ex_gt0{ii}./E_max0{ii};
    Ey_gt{ii} = Ey_gt0{ii}./E_max0{ii};
    Ez_gt{ii} = Ez_gt0{ii}./E_max0{ii};
end
I_gt = I_gt0./max(I_gt0(:));

% retrieved PSF
Ex_ret = cell(1,n_beam);
Ey_ret = cell(1,n_beam);
Ez_ret = cell(1,n_beam);
Ex_delta = cell(1,n_beam);
I_ret_all = cell(1,n_beam);
I_ret = zeros(psfResHi,psfResHi,lz);
error_E_tmp = zeros(1,n_beam);
for ii = 1:n_beam
    [Ex_ret0_tmp,Ey_ret0_tmp,Ez_ret0_tmp] = singleobjectivepsf_noT_MM(Obj,Beam_ret{ii},ScopeHi,Option);
    E_ret0_tmp = cat(3,Ex_ret0_tmp,Ey_ret0_tmp,Ez_ret0_tmp);

    E_gt = cat(3,Ex_gt{ii},Ey_gt{ii},Ez_gt{ii});
    E_ret = complex_correction(E_ret0_tmp,E_gt);

    Ex_ret{ii} = E_ret(:,:,1:lz);
    Ey_ret{ii} = E_ret(:,:,lz+1:2*lz);
    Ez_ret{ii} = E_ret(:,:,2*lz+1:3*lz);
    I_ret_all{ii} = abs(Ex_ret0_tmp).^2 + abs(Ey_ret0_tmp).^2 + abs(Ez_ret0_tmp).^2;
    I_ret = I_ret + I_ret_all{ii};
    I_ret_all{ii} = I_ret_all{ii}./max(I_ret_all{ii}(:));

    Ex_delta{ii} = abs(Ex_gt{ii} - Ex_ret{ii});
end
I_ret = I_ret./max(I_ret(:));

if global_fig_flag == 1
    psf_gain = reshape([1 1 1 1 1],1,1,[]);  % for low NA
    psf_gain_cell = mat2cell(repmat(psf_gain,[1,2+2*n_beam]),1,ones(1,2+2*n_beam),5);
    var_show_I = {I_gt,I_ret};
    for ii = 1:n_beam
        var_show_I{2*ii+1} = I_gt_all{ii};
        var_show_I{2*ii+2} = I_ret_all{ii};
    end
    var_show_I = cellfun(@times,psf_gain_cell,var_show_I,'UniformOutput',false);
    labely = {'|E_{gt}|^2','|E_{ret}|^2'};
    for ii = 1:n_beam
        labely{2*ii+1} = ['|E',num2str(ii),'_{gt}|^2'];
        labely{2*ii+2} = ['|E',num2str(ii),'_{ret}|^2'];
    end
    fig_psfi = figure('Name','psf intensity');
    t_psfi = tiledlayout(fig_psfi,2+2*n_beam,lz,'TileSpacing','compact','Padding','compact');

    for jj = 1:2+2*n_beam
        for ii = 1:lz
            figure(fig_psfi);
            varLz = var_show_I{jj};
            nexttile(t_psfi);
            imagesc(varLz(:,:,ii));axis image xy off;colormap(gca,cmap_PSF);
            clim([0 1]);
            if jj == 1
                title(['z = ',num2str(ScopeHi.zs(ii))]);
            end
            if ii == 1
                axis on;
                xticks([]);
                yticks([]);
                ylabel(labely{jj});
            end
        end
    end

    roi_size = psfResHi/2;
    roi = (psfResHi-roi_size)/2+1:(psfResHi+roi_size)/2;

    % roi = 1:psfResHi;
    idx_z0 = round(lz/2);

    E_gt_roi = cell(1,n_beam);
    E_ret_roi = cell(1,n_beam);
    var_show_I_roi = cell(1,2*n_beam);
    labely = cell(1,2*n_beam);
    for ii = 1:n_beam
        E_gt_roi{ii} = cat(3,Ex_gt{ii}(roi,roi,idx_z0),Ey_gt{ii}(roi,roi,idx_z0),Ez_gt{ii}(roi,roi,idx_z0));
        E_ret_roi{ii} = cat(3,Ex_ret{ii}(roi,roi,idx_z0),Ey_ret{ii}(roi,roi,idx_z0),Ez_ret{ii}(roi,roi,idx_z0));
        var_show_I_roi{2*ii-1} = abs(E_gt_roi{ii}).^2;
        var_show_I_roi{2*ii} = abs(E_ret_roi{ii}).^2;
        labely{2*ii-1} = ['gt ',num2str(ii)];
        labely{2*ii} = ['ret ',num2str(ii)];
    end
    titles = {'|Ex|^2', '|Ey|^2', '|Ez|^2'};

    fig_EFAR = figure('Name','electric field intensity roi');
    t_EFAR = tiledlayout(fig_EFAR,2*n_beam,3,'TileSpacing','compact','Padding','compact');
    for jj = 1:2*n_beam
        for ii = 1:3
            figure(fig_EFAR);
            varLz = var_show_I_roi{jj};
            nexttile(t_EFAR);
            imagesc(varLz(:,:,ii));axis image xy off;colormap(gca,cmap_PSF);
            clim([0 1]);
            if jj == 1
                title(titles{ii});
            end
            if ii == 1
                axis on;
                xticks([]);
                yticks([]);
                ylabel(labely{jj});
            end
        end
    end

    var_show_P_roi = cell(1,2*n_beam);
    for ii = 1:n_beam
        var_show_P_roi{2*ii-1} = angle2pi(E_gt_roi{ii});
        var_show_P_roi{2*ii} = angle2pi(E_ret_roi{ii});
    end
    titles = {'arg(Ex)', 'arg(Ey)', 'arg(Ez)'};

    fig_EFPR = figure('Name','electric field phase roi');
    t_EFPR = tiledlayout(fig_EFPR,2*n_beam,3,'TileSpacing','compact','Padding','compact');
    for jj = 1:2*n_beam
        for ii = 1:3
            figure(fig_EFPR);
            varLz = var_show_P_roi{jj};
            nexttile(t_EFPR);
            imagesc(varLz(:,:,ii));axis image xy off;colormap(gca,cmap_PSFP);
            clim([0 2*pi]);
            if jj == 1
                title(titles{ii});
            end
            if ii == 1
                axis on;
                xticks([]);
                yticks([]);
                ylabel(labely{jj});
            end
        end
    end

    for kk = 1:n_beam
        var_show_E = {abs(Ex_gt{kk}).^2,abs(Ex_ret{kk}).^2,abs(Ey_gt{kk}).^2,abs(Ey_ret{kk}).^2,abs(Ez_gt{kk}).^2,abs(Ez_ret{kk}).^2};
        labely = {'|Ex_{gt}|^2','|Ex_{ret}|^2','|Ey_{gt}|^2','|Ey_{ret}|^2',...
            '|Ez_{gt}|^2','|Ez_{ret}|^2'};

        fig_EFA = figure('Name',['electric field intensity ',num2str(kk)]);
        t_EFA = tiledlayout(fig_EFA,6,lz,'TileSpacing','compact','Padding','compact');
        for jj = 1:6
            for ii = 1:lz
                figure(fig_EFA);
                varLz = var_show_E{jj};
                nexttile(t_EFA);
                imagesc(varLz(:,:,ii));axis image xy off;colormap(gca,cmap_PSF);
                clim([0 1]);
                if jj == 1
                    title(['z = ',num2str(ScopeHi.zs(ii))]);
                end
                if ii == 1
                    axis on;
                    xticks([]);
                    yticks([]);
                    ylabel(labely{jj});
                end
            end
        end

        fig_EFP = figure('Name',['electric field phase ',num2str(kk)]);
        t_EFP = tiledlayout(fig_EFP,6,lz,'TileSpacing','compact','Padding','compact');
        var_show_P = {angle2pi(Ex_gt{kk}),angle2pi(Ex_ret{kk}),angle2pi(Ey_gt{kk}),angle2pi(Ey_ret{kk}),...
            angle2pi(Ez_gt{kk}),angle2pi(Ez_ret{kk})};
        labely = {'arg(Ex_{gt})','arg(Ex_{ret})','arg(Ey_{gt})','arg(Ey_{ret})',...
            'arg(Ez_{gt})','arg(Ez_{ret})'};
        for jj = 1:6
            for ii = 1:lz
                varLz = var_show_P{jj};
                nexttile(t_EFP);
                imagesc(varLz(:,:,ii));axis image xy off;colormap(gca,cmap_PSFP);
                clim([0, 2*pi]);
                if jj == 1
                    title(['z = ',num2str(ScopeHi.zs(ii))]);
                end
                if ii == 1
                    axis on;
                    xticks([]);
                    yticks([]);
                    ylabel(labely{jj});
                end
            end
        end
    end
end
