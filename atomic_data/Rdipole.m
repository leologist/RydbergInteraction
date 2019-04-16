function y=Rdipole(n1,L1,J1,n2,L2,J2)
% Calculates the radial dipole matrix element between states
% |n1 L1 J1> and |n2 L2 J2> in units of a.u.
%
% Vectorized: n1 and n2 should have the same size. The others should wither
% have the same size as n1 (or n2), or be scalar.
%
% For QA purposes, the function usually results in an error message if the pair 
% state does not have an allowed electric dipole transission.
%
% Anlytic formula adapted from [ref 5+7] (with corrections and tweaking)
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.
%

Zcore=1;
if isscalar(L1), L1=n1*0+L1; end; if isscalar(J1), J1=n1*0+J1; end;
if isscalar(L2), L2=n2*0+L2; end; if isscalar(J2), J2=n2*0+J2; end;

nstar1=n1; nstar1(J1>0)=nstar(n1(J1>0),L1(J1>0),J1(J1>0));
nstar2=n2; nstar2(J2>0)=nstar(n2(J2>0),L2(J2>0),J2(J2>0));

I_flip=find(nstar2<nstar1);
temp=    n1(I_flip);    n1(I_flip)=    n2(I_flip);    n2(I_flip)=temp;
temp=    L1(I_flip);    L1(I_flip)=    L2(I_flip);    L2(I_flip)=temp;
temp=nstar1(I_flip);nstar1(I_flip)=nstar2(I_flip);nstar2(I_flip)=temp;

omega=-(1./nstar2.^2-1./nstar1.^2)/2;

dn=nstar2-nstar1;
ncstar=((nstar1.^(-2)+nstar2.^(-2))/2).^(-1/2);
gamma_star=omega.*ncstar.^3;
m=round(gamma_star-dn);
gamma=m+dn;
nc=(gamma./omega).^(1/3);

lambdac=(L1+L2)/2+1/2; eta=lambdac./nc; ec=sqrt(1-eta.^2);

x=gamma.*ec;
sincs=sin(pi*gamma)./(pi*gamma);

dL=L2-L1;
h=m+(n2-n1)-dL;

Dp=(-angerj(gamma,-x,1)+dL.*eta./ec.*(angerj(gamma,-x,0)-sincs))./gamma;
Dr=Dp+(1-ec).*sincs;

y= (-1).^(h+1) .* nc.^5./Zcore./((nstar1.*nstar2).^(3/2)).*Dr;

% References
% [5] Kaulakys, JPhysB, 28, 4963-4971 (1995)
% [7] Kamta,... Oumarou, JPhysB 31 (1998) 963–997.