function [Ex,Ey,Ez] = simple_singleobjectivepsf_MM(sys,PHI,beam)
%SIMPLE_SINGLEOBJECTIVEPSF_CZT is the simplified version of
%'singleobjectivepsf.m'
E0 = beam.amp.*exp(1i*beam.phs);

[Px,Py,Pz] = coortrans(beam.plr(:,:,1),beam.plr(:,:,2),...
    0,sys.theta,sys.phi,'o2i');

E0 = sys.prefix.*E0;
Ex0 = E0.*Px;
Ey0 = E0.*Py;
Ez0 = E0.*Pz;

Ewx = Ex0.*PHI;
Ewy = Ey0.*PHI;
Ewz = Ez0.*PHI;

Ex = pagemtimes(pagemtimes(sys.My,Ewx),sys.Mx);
Ey = pagemtimes(pagemtimes(sys.My,Ewy),sys.Mx);
Ez = pagemtimes(pagemtimes(sys.My,Ewz),sys.Mx);

% energy normalization
Ex = Ex.*sys.df^2;
Ey = Ey.*sys.df^2;
Ez = Ez.*sys.df^2;
end