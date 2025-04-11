function [beam, I0, grad, Option] = vecPrPar_nbeam(sys,beam,Option)
%vecPrPar_nbeam implements parallel phase retrieval (vectorial version)

% forward focusing
Ex0 = cell(1,Option.n_beam);
Ey0 = cell(1,Option.n_beam);
Ez0 = cell(1,Option.n_beam);
Ix0 = zeros(sys{1}.ly,sys{1}.lx,'like',sys{1}.Ir);
Iy0 = Ix0;
Iz0 = Ix0;

for ii = 1:Option.n_beam
    [Ex0{ii},Ey0{ii},Ez0{ii}] = simple_singleobjectivepsf_MM(sys{ii},sys{ii}.PHI,beam{ii});

    Ix0 = Ix0 + abs(Ex0{ii}).^2;
    Iy0 = Iy0 + abs(Ey0{ii}).^2;
    Iz0 = Iz0 + abs(Ez0{ii}).^2;
end
switch Option.Channel_num
    case 1
        I0 = Option.w_x*Ix0 + Option.w_y*Iy0 +Option.w_z*Iz0;
    case 2
        I0_xplr = Option.w_x_xplr*Ix0 + Option.w_y_xplr*Iy0 + Option.w_z_xplr*Iz0;
        I0_yplr = Option.w_x_yplr*Ix0 + Option.w_y_yplr*Iy0 + Option.w_z_yplr*Iz0;
        I0 = cat(3,I0_xplr,I0_yplr);
end

% gradient of loss function w.r.t. intensity
switch Option.Channel_num
    case 1
        pL_pI = 2*(I0 - sys{1}.Target.TI);
    case 2
        pL_pIx = 2*(I0_xplr - sys{1}.Target.TI_xplr);
        pL_pIy = 2*(I0_yplr - sys{1}.Target.TI_yplr);
end

batch_num = size(Ix0,3); lr_w = Option.lr_w/batch_num;
grad = cell(1,Option.n_beam);

for ii = 1:Option.n_beam
    switch Option.Channel_num
        case 1
            pI_pEx = 2*Option.w_x*Ex0{ii};
            pI_pEy = 2*Option.w_y*Ey0{ii};
            pI_pEz = 2*Option.w_z*Ez0{ii};
            
            pL_pEx = pL_pI.*pI_pEx;
            pL_pEy = pL_pI.*pI_pEy;
            pL_pEz = pL_pI.*pI_pEz;
            
            pI_pWx = abs(Ex0{ii}).^2; pI_pWy = abs(Ey0{ii}).^2; pI_pWz = abs(Ez0{ii}).^2;
            
            Option.w_x = Option.w_x - lr_w*sum(pL_pI.*pI_pWx,'all');
            Option.w_y = Option.w_y - lr_w*sum(pL_pI.*pI_pWy,'all');
            Option.w_z = Option.w_z - lr_w*sum(pL_pI.*pI_pWz,'all');
            Option.w_x = max([Option.w_x,0]);
            Option.w_y = max([Option.w_y,0]);
            Option.w_z = max([Option.w_z,0]);
            
        case 2
            pIx_pEx = 2*Option.w_x_xplr*Ex0{ii};
            pIx_pEy = 2*Option.w_y_xplr*Ey0{ii};
            pIx_pEz = 2*Option.w_z_xplr*Ez0{ii};
            pIy_pEx = 2*Option.w_x_yplr*Ex0{ii};
            pIy_pEy = 2*Option.w_y_yplr*Ey0{ii};
            pIy_pEz = 2*Option.w_z_yplr*Ez0{ii};

            pL_pEx = pL_pIx.*pIx_pEx + pL_pIy.*pIy_pEx;
            pL_pEy = pL_pIx.*pIx_pEy + pL_pIy.*pIy_pEy;
            pL_pEz = pL_pIx.*pIx_pEz + pL_pIy.*pIy_pEz;

            pI_pWx = abs(Ex0{ii}).^2;
            pI_pWy = abs(Ey0{ii}).^2;
            pI_pWz = abs(Ez0{ii}).^2;
            
            Option.w_x_xplr = Option.w_x_xplr - lr_w*sum(pL_pIx.*pI_pWx,'all');
            Option.w_y_xplr = Option.w_y_xplr - lr_w*sum(pL_pIx.*pI_pWy,'all');
            Option.w_z_xplr = Option.w_z_xplr - lr_w*sum(pL_pIx.*pI_pWz,'all');
            Option.w_x_yplr = Option.w_x_yplr - lr_w*sum(pL_pIy.*pI_pWx,'all');
            Option.w_y_yplr = Option.w_y_yplr - lr_w*sum(pL_pIy.*pI_pWy,'all');
            Option.w_z_yplr = Option.w_z_yplr - lr_w*sum(pL_pIy.*pI_pWz,'all');

            Option.w_x_xplr = max([Option.w_x_xplr,0]);
            Option.w_y_xplr = max([Option.w_y_xplr,0]);
            Option.w_z_xplr = max([Option.w_z_xplr,0]);
            Option.w_x_yplr = max([Option.w_x_yplr,0]);
            Option.w_y_yplr = max([Option.w_y_yplr,0]);
            Option.w_z_yplr = max([Option.w_z_yplr,0]);
    end
    
    PSFx.amp = abs(pL_pEx); PSFy.amp = abs(pL_pEy); PSFz.amp = abs(pL_pEz);
    PSFx.phs = angle(pL_pEx); PSFy.phs = angle(pL_pEy); PSFz.phs = angle(pL_pEz);
        
    % backward retrieval
    [ampX,phsX] = simple_PR_MM(sys{ii},sys{ii}.PHI,PSFx);
    [ampY,phsY] = simple_PR_MM(sys{ii},sys{ii}.PHI,PSFy);
    [ampZ,phsZ] = simple_PR_MM(sys{ii},sys{ii}.PHI,PSFz);

    Ewx = ampX.*exp(1i*phsX);
    Ewy = ampY.*exp(1i*phsY);
    Ewz = ampZ.*exp(1i*phsZ);

    % inverse projection (from focus to pupil plane)
    [grad_x,grad_y,~] = coortrans(Ewx,Ewy,Ewz,sys{ii}.theta,sys{ii}.phi,'i2o');
    grad{ii} = cat(3,grad_x,grad_y);

    gx = beam{ii}.amp.*exp(1i*beam{ii}.phs).*beam{ii}.plr(:,:,1);
    gy = beam{ii}.amp.*exp(1i*beam{ii}.phs).*beam{ii}.plr(:,:,2);
    Epx = gx - Option.step_size*grad_x;
    Epy = gy - Option.step_size*grad_y;
    
    % remove checkboard pattern on retrieved pupil
    if length(ampX) > 100
        px_eff_real_filt = fftconv3(real(Epx),Option.h);
        px_eff_imag_filt = fftconv3(imag(Epx),Option.h);
        py_eff_real_filt = fftconv3(real(Epy),Option.h);
        py_eff_imag_filt = fftconv3(imag(Epy),Option.h);
        Epx = px_eff_real_filt + 1i*px_eff_imag_filt;
        Epy = py_eff_real_filt + 1i*py_eff_imag_filt;
    end
    
    beam{ii} = pupilregen(beam{ii},Epx,Epy,sys{ii});
end

end
