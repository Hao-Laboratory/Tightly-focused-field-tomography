function C = pagemtimes(A,B)
%PAGEMTIMES is the same as pagemtimes in higher MATLAB version
%   
% LIU Xin
% liuxin2018@zju.edu.cn
% Mar.10, 2024

nrow = size(A,1);
ncol = size(B,2);
np_A = size(A,3);
np_B = size(B,3);
npag = max(np_A,np_B);

C = zeros(nrow,ncol,npag,'like',A);
for ii = 1:npag
    C(:,:,ii) = A(:,:,min(ii,np_A))*B(:,:,min(ii,np_B));
end
end

