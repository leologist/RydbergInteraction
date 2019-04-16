function Vdd = getVdd(N1, L1, J1, M1, N2, L2, J2, M2, geom)
%Vdd = getVdd(N1, L1, J1, M1, N2, L2, J2, M2, geom)
%   Calculate dipole matrix element between two molecular states
%   i.e. <N1 L1 J1 M1|V_dipole|N2 L2 J2 M2>=Vdd/R^3
%   return Vdd in a.u. (2*Ry*a0^3)
%
%   Ni=[n n'] Li=[l l'] Ji=[j j'] Mi=[m m']
%
%   Currently accepts only:
%   geom.type='free space', for which geom.theta is the angle between the dipoles.
% 
%
%   Also possible to run with
%   Vdd = getVdd(N1, L1, J1, M1, geom)
%   where we implicitly compute X2 = fliplr(X1)
%
%   allows any input, and any dipole-forbidden transitions yield Vdd=0
%
%   somewhat vectorized (X1 and/or X2 can be k-by-2)
%
%   author: Leo Zhou

%% parsing inputs
if nargin == 5
    geom = N2;
    N2 = fliplr(N1);
    L2 = fliplr(L1);
    J2 = fliplr(J1);
    M2 = fliplr(M1);
end

if size(N1,1) == 1 & size(N2, 1) > 1
    one1 = ones(size(N2,1), 1);
    N1 = one1*N1; L1 = one1*L1; J1 = one1*J1; M1 = one1*M1;
elseif size(N2,1) == 1 & size(N1, 1) > 1
    one1 = ones(size(N1,1), 1);
    N2 = one1*N2; L2 = one1*L2; J2 = one1*J2; M2 = one1*M2;
end

sizes = diff([size(N1);size(N2);size(L1);size(L2);size(J1);size(J2);size(M1);size(M2)]);
if any(sizes(:))
    error('Array size input invalid');
end

%%

dp = cell(2,1);
for pos = [2,1]
    ind_allowed = find((abs(L2(:,pos) - L1(:,pos)) == 1) &...
        (abs(J2(:,pos)-J1(:,pos)) <= 1) & (abs(M2(:,pos)-M1(:,pos)) <= 1));
    dp{pos} = zeros(size(N2,1),1);
    if ~isempty(ind_allowed)
        dp{pos}(ind_allowed) = ...
            Rdipole_table(N1(ind_allowed,pos),L1(ind_allowed,pos),J1(ind_allowed,pos),...
            N2(ind_allowed,pos), L2(ind_allowed,pos), J2(ind_allowed,pos))...
            .* reduced_LJM(L1(ind_allowed,pos),J1(ind_allowed,pos),M1(ind_allowed,pos),...
            L2(ind_allowed,pos), J2(ind_allowed,pos), M2(ind_allowed,pos));
    end
end
dp = (dp{1}.*dp{2});

dm_L = M2(:,1)-M1(:,1);
dm_R = M2(:,2)-M1(:,2);
I1p=(dm_L==1);I1z=(dm_L==0);I1m=(dm_L==-1);
I2p=(dm_R==1);I2z=(dm_R==0);I2m=(dm_R==-1);

theta=geom.angle;
geom_factor_free_space = I1p.*I2m + I1m.*I2p + (I1z.*I2z)*(1-3*cos(theta)^2)-...
    (3/2)*sin(theta)^2*(I1p.*I2p+I1p.*I2m+I1m.*I2p+I1m.*I2m)-...
    (3/sqrt(2))*sin(theta)*cos(theta)*(I1p.*I2z+I1m.*I2z+I1z.*I2p+I1z.*I2m); % ref. [9]
Vdd = geom_factor_free_space.*dp;


end

%% References
% [9] Reinhard, ..., Raithel PRA 75, 032712 (2007).