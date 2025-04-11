function [x_cent,y_cent,x_size,y_size] = ...
    psf_cent_size(psfDenoised,th,pixelSize,method_flag)
%PSF_CENT_SIZE calculates the specification of PSF, including size and
%center
%   
% Xin Liu
% liuxin2018@zju.edu.cn
% Nov. 22, 2023

sz = size(psfDenoised,[1,2]);
Sx = pixelSize*(sz(2)-1); Sy = pixelSize*(sz(1)-1);
x = linspace(-Sx/2,Sx/2,sz(2));
y = linspace(-Sy/2,Sy/2,sz(1));

psfDenoised = psfDenoised./max(psfDenoised(:));

if nargin == 3
    method_flag = 'centriod';
end

switch method_flag
    case 'intuitive'  % preferred for computational imaging since this determine the psf size according to the noise floor
        % th = 1e-2;  % preferred value
        S_row = 1:sz(1);
        S_col = 1:sz(2);
        [cols,rows] = meshgrid(S_col,S_row);

        % find centriod
        centCol = round(sum(cols.*psfDenoised,"all")./sum(psfDenoised,"all"));
        centRow = round(sum(rows.*psfDenoised,"all")./sum(psfDenoised,"all"));

        mask_temp = double(psfDenoised>th);
        I_temp = imfill(mask_temp,'holes');  % fill holes

        [~,max_idx] = max(psfDenoised(:));  % obtain the index of the peak
        I_temp = grayconnected(I_temp,centRow,centCol);  % connect discrete regions
        I_temp = double(I_temp);
        if I_temp(max_idx) ~= 1
            I_temp = 1 - I_temp;
        end
        mask = grayconnected(I_temp,centRow,centCol);  % still may be failed if 'th' is too large
        % figure;
        % imagesc(x,y,mask);axis image xy;colorbar;

        [row_idx,col_idx] = find(mask);
        x_list = -Sx/2 + (col_idx-1).*pixelSize;
        y_list = -Sy/2 + (row_idx-1).*pixelSize;
        x_min = min(x_list); x_max = max(x_list);
        y_min = min(y_list); y_max = max(y_list);

        x_cent = 0.5*(x_min+x_max);
        y_cent = 0.5*(y_min+y_max);

    case 'centriod'  % centroid method
        % th = 1e-2;  % preferred value
        [xx,yy] = meshgrid(x,y);
        x_cent = sum(xx.*psfDenoised,"all")./sum(psfDenoised,"all");
        y_cent = sum(yy.*psfDenoised,"all")./sum(psfDenoised,"all");

        energy_col = cumsum(sum(psfDenoised,1))./sum(psfDenoised,"all");
        energy_row = cumsum(sum(psfDenoised,2))./sum(psfDenoised,"all");

        idx_x_min = find(energy_col>th,1,'first');
        idx_x_max = find(energy_col>1-th,1,'first');
        idx_y_min = find(energy_row>th,1,'first');
        idx_y_max = find(energy_row>1-th,1,'first');
        x_min = x(idx_x_min); x_max = x(idx_x_max);
        y_min = y(idx_y_min); y_max = y(idx_y_max);
    case 'rms'  % root-mean-square method
        sf = 20;  % scale factor
        [xx,yy] = meshgrid(x,y);
        x_cent = sum(xx.*psfDenoised,"all")./sum(psfDenoised,"all");
        y_cent = sum(yy.*psfDenoised,"all")./sum(psfDenoised,"all");
        sigma_x_sq = sf*sum((xx-x_cent).^2.*psfDenoised,"all")./sum(psfDenoised,"all");
        sigma_y_sq = sf*sum((yy-y_cent).^2.*psfDenoised,"all")./sum(psfDenoised,"all");
        x_max = x_cent+0.5*sqrt(sigma_x_sq); x_min = x_cent-0.5*sqrt(sigma_x_sq);
        y_max = y_cent+0.5*sqrt(sigma_y_sq); y_min = y_cent-0.5*sqrt(sigma_y_sq);
end

% figure;
% imagesc(x,y,psfDenoised.^1);axis image xy;colorbar;
% hold on;
% plot([x_min,x_max,x_max,x_min,x_min],[y_min,y_min,y_max,y_max,y_min],'LineWidth',1);

x_size = x_max - x_min;
y_size = y_max - y_min;
end

