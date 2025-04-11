function plrOut = plrBasisTransform(plrIn,flag)
%PLRBASISTRANSFORM convert polarization state between two bases
%   
% LIU Xin
% liuxin2018@zju.edu.cn
% May. 13, 2022

pupilRes = size(plrIn,1);
xy = linspace(-1,1,pupilRes);
[xx,yy] = meshgrid(xy,xy);
[phi,~] = cart2pol(xx,yy);

plrOut = zeros(size(plrIn));
switch flag
    case 'xy2circle'  
        % P_x = plrIn(:,:,1);  P_y = plrIn(:,:,2);
        % plrOut(:,:,1) = P_left;  plrOut(:,:,2) = P_right;
        plrOut(:,:,1) = 1/sqrt(2)*(plrIn(:,:,1) - 1i*plrIn(:,:,2));
        plrOut(:,:,2) = 1/sqrt(2)*(plrIn(:,:,1) + 1i*plrIn(:,:,2));
    case 'xy2cylinder'  
        % P_x = plrIn(:,:,1);  P_y = plrIn(:,:,2);
        % plrOut(:,:,1) = P_rho;  plrOut(:,:,2) = P_phi;
        plrOut(:,:,1) =  plrIn(:,:,1).*cos(phi) + plrIn(:,:,2).*sin(phi);
        plrOut(:,:,2) = -plrIn(:,:,1).*sin(phi) + plrIn(:,:,2).*cos(phi);
    case 'circle2xy'
        % P_left = plrIn(:,:,1);  P_right = plrIn(:,:,2);
        % plrOut(:,:,1) = P_x;  plrOut(:,:,2) = P_y;
        plrOut(:,:,1) = 1 /sqrt(2)*(plrIn(:,:,1) + plrIn(:,:,2));
        plrOut(:,:,2) = 1i/sqrt(2)*(plrIn(:,:,1) - plrIn(:,:,2));
    case 'circle2cylinder'  
        % P_left = plrIn(:,:,1);  P_right = plrIn(:,:,2);
        % plrOut(:,:,1) = P_rho;  plrOut(:,:,2) = P_phi;
        plrOut(:,:,1) = 1 /sqrt(2)*(plrIn(:,:,1).*exp(1i*phi) + plrIn(:,:,2).*exp(-1i*phi));
        plrOut(:,:,2) = 1i/sqrt(2)*(plrIn(:,:,1).*exp(1i*phi) - plrIn(:,:,2).*exp(-1i*phi));
    case 'cylinder2xy'  
        % P_rho = plrIn(:,:,1);  P_phi = plrIn(:,:,2);
        % plrOut(:,:,1) = P_x;  plrOut(:,:,2) = P_y;
        plrOut(:,:,1) = plrIn(:,:,1).*cos(phi) - plrIn(:,:,2).*sin(phi);
        plrOut(:,:,2) = plrIn(:,:,1).*sin(phi) + plrIn(:,:,2).*cos(phi);
    case 'cylinder2circle'  
        % P_rho = plrIn(:,:,1);  P_phi = plrIn(:,:,2);
        % plrOut(:,:,1) = P_left;  plrOut(:,:,2) = P_right;
        plrOut(:,:,1) = 1/sqrt(2)*(plrIn(:,:,1) - 1i*plrIn(:,:,2)).*exp(-1i*phi);
        plrOut(:,:,2) = 1/sqrt(2)*(plrIn(:,:,1) + 1i*plrIn(:,:,2)).*exp( 1i*phi);
    otherwise
        error('No such an option, please try again.');
end
end

