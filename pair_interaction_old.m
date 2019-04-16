function [y, energies]=pair_interaction_old(N1,L1,J1,M1,geom,R)
% Calculates the dipole-dipole and VdV interaction of the pair state:
% N1=[n n'] L1=[l l'] J1=[j j'] M1=[m m'],
% by going over many relevant "virtual" pair states (N2s,L2s,J2s,M2s).
% 
% R is in units of [a0] (a.u.)
% if R is not given, the function caluclates C6
%
% Currently accepts only:
% geom.type='free space', for which geom.theta is the angle between the dipoles.
%
% Not vectorized.
% Note: might "crush" the matlab if a pair state is given as an input, for which there is 
% an electric-dipole allowed transition.
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.
%
% Improved by Leo Zhou, Harvard University, 9/16/2016

global ReducedMassFactor % =(1+me/ma)
if isempty(ReducedMassFactor), ReducedMassFactor=1; end;

if ~exist('N2s_distance','var')
    N2s_distance=7; % determines how far in the n-manifold to consider states for the virtual pair.
end
threshold = 1e-6;

N2=cell(2,1);M2=N2;dE=N2;dp=N2;
for r=[1 2]
    % going over all relevant allowed "virtual" states
    % JDT: construct arrays of all dipole allowed virtual states, without
    % regard for whether these states are physical. Change L by one, have
    % J' = L'+/- 1/2,
    % and M by +-1 and 0.
    N2s=N1(r)+(-N2s_distance:N2s_distance);
    L2s=L1(r)  +[-[1 1 1   1 1 1]  [1 1 1   1 1 1]];
    J2s=L2s+0.5*[-[1 1 1] [1 1 1] -[1 1 1] [1 1 1]];
    M2s=M1(r)  +[ -1 0 1  -1 0 1   -1 0 1  -1 0 1 ];
    % The "forbidden" are: (i) dJ=0,+1,-1   (ii) J2 not negative, and can support M2 (iii) L2 not negative
    % JDT: Now remove the states that are not physical (in terms of angular
    % momentum quantum numbers only)
    Iforbid=(abs(J2s-J1(r))>1 | J2s<abs(M2s) | L2s<0);L2s(Iforbid)=[];J2s(Iforbid)=[];M2s(Iforbid)=[];
    
    % JDT: now construct arrays of angular momentum quantum numbers for
    % each principal quantum number
    [kuku,L2s]=meshgrid(N2s,L2s);[kuku,J2s]=meshgrid(N2s,J2s);[N2s,M2s]=meshgrid(N2s,M2s); %#ok<ASGLU>
    N2s=N2s(:);J2s=J2s(:);L2s=L2s(:);M2s=M2s(:); % JDT: Flatten into one-dim array
    
    % JDT: now remove unphysical states where L>N
    IlargeL=(L2s>=N2s);N2s(IlargeL)=[];L2s(IlargeL)=[];J2s(IlargeL)=[];M2s(IlargeL)=[];
    
    % JDT: get angular part of matrix element
    reduced_angular=reduced_LJM(L1(r),J1(r),M1(r),L2s,J2s,M2s);
    
    % JDT: get energy difference in AU
    dE{r}=(nstar(N1(r),L1(r),J1(r)).^(-2)-nstar(N2s,L2s,J2s).^(-2))/2/ReducedMassFactor; % To convert to a.u. (factor 2 because the Hartree energy is 2*Ry);
    N1s=N2s*0+N1(r);
    Rdipole12=Rdipole_table(N1s,L1(r),J1(r),N2s,L2s,J2s);
    
    % JDT: calculate total matrix elements, and save magnetic sublevels
    dp{r}=abs(reduced_angular.*Rdipole12);
    M2{r}=M2s;
end
clear N2s L2s J2s M2s
% Calculates energy defects and dipole-dipole interaction for each of states

% JDT: make table of all possible intermediate state pairs, with energy
% difference and products of matrix elements
[dE{1} dE{2}]=meshgrid(dE{1},dE{2});dE=dE{1}+dE{2};
[dp{1} dp{2}]=meshgrid(dp{1},dp{2});dp=dp{2}.*dp{1};

% JDT: also find the polarization of each transition. dM1 is change in
% magnetic sublevel for first atom, and dM2 for second
[dM1 dM2]=meshgrid(M2{1}-M1(1),M2{2}-M1(2));
I1p=(dM1==1);I1z=(dM1==0);I1m=(dM1==-1);
I2p=(dM2==1);I2z=(dM2==0);I2m=(dM2==-1);

theta=geom.angle;
geom_factor_free_space=I1p.*I2m+I1m.*I2p+(I1z.*I2z)*(1-3*cos(theta)^2)-...
    (3/2)*sin(theta)^2*(I1p.*I2p+I1p.*I2m+I1m.*I2p+I1m.*I2m)-...
    (3/sqrt(2))*sin(theta)*cos(theta)*(I1p.*I2z+I1m.*I2z+I1z.*I2p+I1z.*I2m); % ref. [9]

Ires = find(dE==0);
Ioff = find(dE~=0);
if isempty(R)
    % calculates C6 in a.u.
    y = -sum(sum(abs(geom_factor_free_space.*dp).^2./dE));
    energies = [];
else
    % calculates the pair interaction of a given distance R
    Vdd=geom_factor_free_space.*dp./R^3;
%     dW=(dE-sign(dE).*sqrt(dE.^2+4*abs(Vdd).^2))/2; % old
%     y=sum(sum(dW)); % old
    dW(Ioff) = (dE(Ioff) - sign(dE(Ioff)).*sqrt(dE(Ioff).^2+4*abs(Vdd(Ioff)).^2))/2;
    dW(Ires) = -abs(Vdd(Ires));
    y = sum(dW(:));
    
    if 1 
        % diagonalization of a matrix containing only pair interactions of
        % the original pair with other allowed pair:
        %
        %    [    0   V1(R)  V2(R)  V3(R) ... ]
        %    [ V1*(R)  dE1     0      0    0  ]
        % M= [ V2*(R)   0     dE2     0    0  ]
        %    [ V3*(R)   0      0     dE3   0  ]
        %    [   ...    0      0      0   ... ]
        %
        % In princple, all the off-diagonal zeros above should be accupied
        % with the interaction between the other states, and I decided to
        % take that as zero for now. This actually might result in a large
        % error (say, up to a factor of 2) in the case of resonant
        % dipole-dipole interaction, so should consider filling-up the
        % whole matrix in the future.
        
        I=(abs(dW)>=abs(y)*threshold);
        %I=(abs(dW)>=0);
        dEvec=dE(I); Vddvec=Vdd(I);
        M=diag([0, dEvec]);
        M(1,2:end)=Vddvec';
        M(2:end,1)=Vddvec;
        
        if length(M) > 5e3
            error(sprintf('Matrix too large: dim = %i',length(M)));
        end
        
        [V,D]=eig(M);
        [weight_of_original_pair_state,I]=max(abs(V(1,:))); %#ok<ASGLU>
        y=abs(D(I,I));
        energies = diag(D);
    end
end

%% References
% [9] Reinhard, ..., Raithel PRA 75, 032712 (2007).