function y=reduced_LJM(L1,J1,M1,L2,J2,M2)
%
% Calculates the reduced matrix elements <L1 J1 M1|mu_q|L2 J2 M2>
%
% Vectorized.
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.

one1=ones(max([size(L1);size(J1);size(M1);size(L2);size(J2);size(M2)],[],1));
zer0=one1*0;
if isscalar(L1), L1=L1*one1; end
if isscalar(J1), J1=J1*one1; end
if isscalar(M1), M1=M1*one1; end
if isscalar(L2), L2=L2*one1; end
if isscalar(J2), J2=J2*one1; end
if isscalar(M2), M2=M2*one1; end


S=1/2*one1;
reduced_LL=            (-1).^(L1)       .*wigner3j(L1,zer0, one1,zer0,L2,zer0).*sqrt((2*L1+1).*(2*L2+1));
reduced_JJ=reduced_LL.*(-1).^(L1+J2+S+1).*wigner6j(J1,one1,   J2,  L2,S   ,L1).*sqrt((2*J1+1).*(2*J2+1));

if ~isempty(M2)
    q_pol=M1-M2;
    reduced_MM=reduced_JJ.*(-1).^(L1-M1)    .*wigner3j(J1,-M1, one1,q_pol,J2,M2);
else
    reduced_MM_sqr=reduced_JJ.*0;
    M2=M1-1;I=abs(J2*0+M2)<=J2;q_pol=M1-M2;
    if sum(I)>0, reduced_MM_sqr(I)=reduced_MM_sqr(I)+abs(w3j(J1(I),-M1(I), 1,q_pol(I),J2(I),M2(I))).^2;end
    M2=M1  ;I=abs(J2*0+M2)<=J2;q_pol=M1-M2;
    if sum(I)>0, reduced_MM_sqr(I)=reduced_MM_sqr(I)+abs(w3j(J1(I),-M1(I), 1,q_pol(I),J2(I),M2(I))).^2;end
    M2=M1+1;I=abs(J2*0+M2)<=J2;q_pol=M1-M2;
    if sum(I)>0, reduced_MM_sqr(I)=reduced_MM_sqr(I)+abs(w3j(J1(I),-M1(I), 1,q_pol(I),J2(I),M2(I))).^2;end
    reduced_MM=sqrt(reduced_MM_sqr).*reduced_JJ;
end
y=reduced_MM;
