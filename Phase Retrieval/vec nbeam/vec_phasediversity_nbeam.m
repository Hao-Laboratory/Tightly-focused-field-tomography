function [amp, phs, plr, Ir, metrics_val, sys, Option] = vec_phasediversity_nbeam(beam,obj,scope,...
    Target,PhsDiv,Threshold,Option)
% main function for vectorial tomography
%
% *************************************************************************
% Xin Liu
% liuxin.optics@gmail.com

%% built-in value
n_beam = length(beam);
if Option.PltFig == 1
    fig_ite = figure('Name','Vectorial Phase Retrieval','Units','normalized');
    Pos = get(fig_ite,'Position');
    set(fig_ite,'Position',[Pos(1:2),0.5,Pos(4)]);
    t_ite = tiledlayout(fig_ite,n_beam+1,3,'TileSpacing','compact','Padding','none');
end

%% parameters register
sys = cell(1,n_beam);
PHI_temp = cell(1,n_beam);

for ii = 1:n_beam
    [sys{ii}, beam{ii}] = SysInit(obj,beam{ii},scope,Target,Option);
    sys{ii}.prMode = Option.PrMode;
    PHI_temp{ii} = exp(1i*PhsDiv{ii});
end

%% initial guess of pupil function
Option.reini_count = 0;
for ii = 1:n_beam
    beam{ii} = pupilreini(beam{ii},Option);
end

%% initial guess of weight for dipole emission
Option.n_beam = n_beam;
if Option.lr_w~=0
    switch Option.Channel_num
        case 1
            Option.w_x = Option.w_x + 0.05;
            Option.w_y = Option.w_y + 0.05;
            Option.w_z = Option.w_z + 0.05;
        case 2
            Option.w_x_xplr = Option.w_x_xplr + 0.05;
            Option.w_y_xplr = Option.w_y_xplr - 0.05;
            Option.w_z_xplr = Option.w_z_xplr + 0.05;
            
            Option.w_x_yplr = Option.w_x_yplr + 0.05;
            Option.w_y_yplr = Option.w_y_yplr - 0.05;
            Option.w_z_yplr = Option.w_z_yplr + 0.05;
    end
end

switch Option.Channel_num
    case 1
        w_total_1ch = Option.w_x + Option.w_y + Option.w_z;
        Option.w_x = Option.w_x/w_total_1ch;
        Option.w_y = Option.w_y/w_total_1ch;
        Option.w_z = Option.w_z/w_total_1ch;
        
        w_x = ones(1,1);
        if Option.UseGpu == 1
            w_x = gpuArray(w_x);
        end
        w_y = w_x;
        w_z = w_x;
        
    case 2
        w_total_2ch = Option.w_x_xplr + Option.w_y_xplr + Option.w_z_xplr + ...
            Option.w_x_yplr + Option.w_y_yplr + Option.w_z_yplr;
        
        Option.w_x_xplr = Option.w_x_xplr/w_total_2ch;
        Option.w_y_xplr = Option.w_y_xplr/w_total_2ch;
        Option.w_z_xplr = Option.w_z_xplr/w_total_2ch;
        
        Option.w_x_yplr = Option.w_x_yplr/w_total_2ch;
        Option.w_y_yplr = Option.w_y_yplr/w_total_2ch;
        Option.w_z_yplr = Option.w_z_yplr/w_total_2ch;
        
        w_x_yplr = ones(1,1);
        if Option.UseGpu == 1
            w_x_yplr = gpuArray(w_x_yplr);
        end
        w_y_yplr = w_x_yplr;
        w_z_yplr = w_x_yplr;
end

%% generate smooth filter
Option.h = fspecial('gaussian',5,0.5);  % smooth potential artifacts induced by PSF denoising
if Option.UseGpu == 1
    Option.h = gpuArray(Option.h);
end

%% prepare iteration
metrics_val = ones(1,1);
pupil_var = ones(1,1);
pupil_var_temp = zeros(1,n_beam);
if Option.UseGpu == 1
    metrics_val = gpuArray(metrics_val);
    pupil_var = gpuArray(pupil_var);
    pupil_var_temp = gpuArray(pupil_var_temp);
end

rand_N = size(PhsDiv{1},3);
switch Option.Channel_num
    case 1
        Target_temp.TI = sys{1}.Target.TI;
    case 2
        Target_temp.TI_xplr = sys{1}.Target.TI_xplr;
        Target_temp.TI_yplr = sys{1}.Target.TI_yplr;
end

%% start iteration
flag = true;
for nstep = 1:Threshold.MaxIter
    if nstep>Threshold.MaxIter-100
        flag = false;
    end
    
    if flag
        batch_idx = randperm(rand_N,Option.PrBatchSize);  % mini-batch SGD
    else
        batch_idx = 1:rand_N;  % GD
    end
    
    for ii = 1:n_beam
        sys{ii}.PHI = PHI_temp{ii}(:,:,batch_idx);
    end
    
    switch Option.Channel_num
        case 1  % direct detection
            sys{1}.Target.TI = Target_temp.TI(:,:,batch_idx);
            I_target = sys{1}.Target.TI;
            
            w_total_1ch = Option.w_x + Option.w_y + Option.w_z;
            
            w_x(nstep) = Option.w_x/w_total_1ch;
            w_y(nstep) = Option.w_y/w_total_1ch;
            w_z(nstep) = Option.w_z/w_total_1ch;
        case 2  % polarization split
            sys{1}.Target.TI_xplr = Target_temp.TI_xplr(:,:,batch_idx);
            sys{1}.Target.TI_yplr = Target_temp.TI_yplr(:,:,batch_idx);
            I_target = cat(3,sys{1}.Target.TI_xplr,sys{1}.Target.TI_yplr);
            
            w_total_2ch = Option.w_x_xplr + Option.w_y_xplr + Option.w_z_xplr + ...
                Option.w_x_yplr + Option.w_y_yplr + Option.w_z_yplr;
            
            w_x_yplr(nstep) = Option.w_x_yplr/w_total_2ch;
            w_y_yplr(nstep) = Option.w_y_yplr/w_total_2ch;
            w_z_yplr(nstep) = Option.w_z_yplr/w_total_2ch;
    end
    
    for ii = 1:n_beam
        sys{ii}.prMode = Option.PrMode;
    end
    
    % gradient computation and update
    [beam,Ir,grad,Option] = vecPrPar_nbeam(sys,beam,Option);

    for jj = 1:n_beam
        grad_sum = sqrt(sum(abs(grad{jj}).^2,3));
        pupil_mask = beam{jj}.amp~=0;
        pupil_var_temp(jj) = mean(grad_sum(pupil_mask))/mean(beam{jj}.amp(pupil_mask));
    end
    pupil_var(nstep) = mean(pupil_var_temp);
    metrics_val(nstep) = 1-rmse(Ir,I_target);

    pupil_var_show = pupil_var;
    metrics_val_show = metrics_val;
    
    % figure
    if Option.PltFig == 1
        if nstep == 1 || mod(nstep,1e2) == 0 || nstep == Threshold.MaxIter
            
            for ii = 1:n_beam
                if contains(Option.PrMode,'amp')
                    nexttile(t_ite,1+(ii-1)*3);
                    pupilshow(beam{ii}.amp);
                    colormap(gca,'gray');
                    colorbar;
                    title('amp');
                end
                
                if contains(Option.PrMode,'phs')
                    nexttile(t_ite,2+(ii-1)*3);
                    pupilshow(wrapTo2Pi(beam{ii}.phs));
                    colormap(gca,'gray');
                    colorbar;
                    caxis([0 2*pi]);
                    title('phs');
                end
                
                % if contains(Option.PrMode,'plr')
                %     nexttile(t_ite,3+(ii-1)*3);
                %     pupilshowplr([],gather(beam{ii}.plr.*beam{ii}.amp),'num',Option.PlrNum,'style',...
                %         Option.PlrStyle,'bgcolor',Option.PlrBgcolor,'cmap',Option.PlrCmap,'normflag',Option.PlrNorm);
                %     colorbar;
                %     title('plr');
                % end
            end
            
            idx_show = randi(size(I_target,3));
            nexttile(t_ite,n_beam*3+1);
            imagesc(I_target(:,:,idx_show));
            axis image xy off;
            colorbar;
            title('gt PSF');
            
            nexttile(t_ite,n_beam*3+2);
            imagesc(Ir(:,:,idx_show));
            axis image xy off;
            colorbar;
            title('retrieved PSF');
            
            nexttile(t_ite,n_beam*3+3);
            plot(1:nstep,metrics_val_show,...
                'DisplayName',['Fidelity: ', sprintf('%.3e',metrics_val_show(nstep))],'LineWidth',1);
            hold on;
            plot(1:nstep,pupil_var_show,...
                'DisplayName',['grad: ', sprintf('%.3e',pupil_var_show(nstep))],'LineWidth',1);
            switch Option.Channel_num
                case 1
                    plot(1:nstep,w_x,'DisplayName', ['w_x: ', sprintf('%.2f',w_x(end))],'LineWidth',1);
                    plot(1:nstep,w_y,'DisplayName', ['w_y: ', sprintf('%.2f',w_y(end))],'LineWidth',1);
                    plot(1:nstep,w_z,'DisplayName', ['w_z: ', sprintf('%.2f',w_z(end))],'LineWidth',1);
                case 2
                    plot(1:nstep,w_x_yplr,'DisplayName',['w_x yplr: ', sprintf('%.2f',w_x_yplr(end))],'LineWidth',1);
                    plot(1:nstep,w_y_yplr,'DisplayName',['w_y yplr: ', sprintf('%.2f',w_y_yplr(end))],'LineWidth',1);
                    plot(1:nstep,w_z_yplr,'DisplayName',['w_z yplr: ', sprintf('%.2f',w_z_yplr(end))],'LineWidth',1);
            end
            hold off;
            xlim([1 Threshold.MaxIter]);
            xticklabels(cellfun(@num2str,num2cell(xticks*1e-3),'UniformOutput',false));
            ylim([0 1]);
            xlabel('Step (k)');
            ylabel('Value');
            grid on;
            legend('NumColumns',2);
            title([num2str(nstep),' iterations']);
            drawnow;
        end
    end
end
fprintf('Converged, phase retrieval finished!.\n');

%% prepare output
ampTotal = [];
for ii = 1:n_beam
    ampTotal = cat(3,ampTotal,beam{ii}.amp);
end

for ii = 1:n_beam
    beam{ii}.amp = beam{ii}.amp./max(ampTotal(:));
end

amp = cell(1,n_beam);
phs = cell(1,n_beam);
plr = cell(1,n_beam);
for ii = 1:n_beam
    beam{ii}.amp = circularpupil(beam{ii}.amp);
    beam{ii}.phs = circularpupil(beam{ii}.phs);
    beam{ii}.plr = circularpupil(beam{ii}.plr);
    amp{ii} = beam{ii}.amp;
    phs{ii} = wrapTo2Pi(beam{ii}.phs);
    plr{ii} = beam{ii}.plr;
    
    if Option.UseGpu == 1
        amp{ii} = gather(amp{ii});
        phs{ii} = gather(phs{ii});
        plr{ii} = gather(plr{ii});
    end
    
    amp{ii} = double(amp{ii});
    phs{ii} = double(phs{ii});
    plr{ii} = double(plr{ii});
end

if Option.UseGpu == 1
    Ir = gather(Ir);
end
Ir = double(Ir);

switch Option.Channel_num
    case 1
        w_total_1ch = Option.w_x + Option.w_y + Option.w_z;
        Option.w_x = Option.w_x/w_total_1ch;
        Option.w_y = Option.w_y/w_total_1ch;
        Option.w_z = Option.w_z/w_total_1ch;
        
    case 2
        w_total_2ch = Option.w_x_xplr + Option.w_y_xplr + Option.w_z_xplr + ...
            Option.w_x_yplr + Option.w_y_yplr + Option.w_z_yplr;
        
        Option.w_x_xplr = Option.w_x_xplr/w_total_2ch;
        Option.w_y_xplr = Option.w_y_xplr/w_total_2ch;
        Option.w_z_xplr = Option.w_z_xplr/w_total_2ch;
        
        Option.w_x_yplr = Option.w_x_yplr/w_total_2ch;
        Option.w_y_yplr = Option.w_y_yplr/w_total_2ch;
        Option.w_z_yplr = Option.w_z_yplr/w_total_2ch;
end
end