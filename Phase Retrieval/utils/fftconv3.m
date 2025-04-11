function convResult=fftconv3(a,b)
% matrix expansion
[a1,a2,a3] = size(a);
[b1,b2,b3] = size(b);

if any([mod(b1,2)==0,mod(b2,2)==0,mod(b3,2)==0])
    warning('Kernel with odd pixel number is recommended.')
end

c = zeros(a1+b1-1,a2+b2-1,a3+b3-1,'like',b);
d = c;

c(1:a1,1:a2,1:a3) = a;
C = fftn(c);
clear c;
d(1:b1,1:b2,1:b3) = b;
D = fftn(d);
clear d;

E = C.*D;
clear C D;

E = ifftn(E);

convResult = E((floor(b1/2)+1):(floor(b1/2)+a1), ...
    (floor(b2/2)+1):(floor(b2/2)+a2), ...
    (floor(b3/2)+1):(floor(b3/2)+a3));
end