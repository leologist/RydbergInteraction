function [y, energies, overlaps, info] = pair_interaction(N1,L1,J1,M1,geom,R,N2s_distance,threshold)
% Calculates the dipole-dipole and VdV interaction of the pair state:
% N1=[n n'] L1=[l l'] J1=[j j'] M1=[m m'],
% by going over many relevant "virtual" pair states (N2s,L2s,J2s,M2s),
% building Hamiltonian and diagonalizing it.
%
% Usage:
% [y, energies, overlaps, info] = pair_interaction(N1,L1,J1,M1,geom,R);
% [y, energies, overlaps, info] = pair_interaction(N1,L1,J1,M1,geom,R,N2s_distance,threshold)
% R is in units of [a0] (a.u.) (R = (true R in meters)/a0)
% Optional parameters: N2_distance (default = 7), and threshold (default = 1e-4)
% Currently accepts only:
% geom.type='free space', for which geom.theta is the angle between the dipoles.
% 
% y = energy shift from a state adiabatically connected to asymptotic state
% energies = eigenvalues from exact diagonalization
% overlaps = eigenstate overlap with asymptotic state (e.g. <sp+ps|E_i>)
% info.overlap_vectors = asymptotic state vector with which overlaps are computed
% info.state_labels = state labels of each index for overlap_vectors
%                     (n1, n2, l1, l2, j1, j2, m1, m2)
% info.HamResSect = truncated Hamiltonian for the resonant sector with R=1;
%
% NOT VECTORIZED
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.
%
% Improved by Leo Zhou, Harvard University, 6 Oct 2016
%   - included states that couples to the states resonant with input state
%       - more accurate for computing flip-flop interaction for |nSnP>
%   - additionally outputs overlaps of eigenstates with asymptotic states
%   - the asymptotic states are eigenstates of info.HamResSect
%   - implemented hashing of state labels via Chinese Remainder Theorem

global ReducedMassFactor % =(1+me/ma)
if isempty(ReducedMassFactor)
    ReducedMassFactor=1;
end

%% program parameters

% determines how far in the n-manifold to consider states for the virtual pair.
if ~exist('N2s_distance','var')
    N2s_distance=7;
end

% only keep states who contributes energy shift larger than threshold*(estimated energy shfit)
if ~exist('threshold', 'var')
    threshold = 1e-4;
end

hashprimes = [157, 151, 7, 11, 13, 17, 19, 23];
hashmultipliers = [25826717917, 57208448297, -75555180987, -80134282865, -40683558993, 72592232713, -74229651496, -61320146888]';
% one-to-one hashing from state labels to single integer
% hash function based on Chinese Remainder Theorem
% hashmultipliers = getCRHashMultipliers(hashprimes);

%% Get all states that couple to 1) input states AND 2) states that resonantly couples to input state
% as well as calculating their Forster energy defect and dipole matrix elements

dEvecs = {}; % a cell for each nontrivial resonant state
Vddvecs = {};
state_hashes = {};

% start with the inputted state
state_hash = getHashInd([N1, L1, J1, M1]); % this is our master list of states
resonantSector = state_hash; % master list of resonant states
dim_ResonantSector = 1; % may discover other states on resonance later

ind = 1;
while ind <= dim_ResonantSector
    state = getStateLabel(resonantSector(ind));
    
    % get Forster defects and dipole matrix elements
    [dEvecs{ind}, Vddvecs{ind}, state_hashes{ind}] = ...
        getInteractions(state(1:2), state(3:4), state(5:6), state(7:8));
    
    Iresonant = (dEvecs{ind} == 0 & Vddvecs{ind} ~=0); % find resonant states
    
    dW=(dEvecs{ind}-(sign(dEvecs{ind})+Iresonant).*sqrt(dEvecs{ind}.^2+4*abs(Vddvecs{ind}).^2))/2;
    
    y = sum(sum(dW));
    
    % keep only states who contribute above threshold, or resonant states
    I = abs(dW)> abs(y)*threshold | Iresonant;
    
    if ~nnz(I) % if no such states
        [~, I] = max(abs(dW(:))); % keep the largest interacting one
    end
    
    dEvecs{ind} = dEvecs{ind}(I); Vddvecs{ind} = Vddvecs{ind}(I);
    state_hashes{ind} = state_hashes{ind}(I);

    % check for new states to add to resonant sector
    Iresonant = (dEvecs{ind} == 0 & Vddvecs{ind} ~=0);
    newResonantStates = state_hashes{ind}(Iresonant);
    [~, Inewres] = setdiff(newResonantStates, resonantSector);
    if ~isempty(Inewres) % add new states to resonant sector
        resonantSector = [resonantSector; newResonantStates(Inewres)];
        dim_ResonantSector = length(resonantSector);
    end
    
    % update master list of states
    [~, Inew] = setdiff(state_hashes{ind}, state_hash);
    state_hash = [state_hash; state_hashes{ind}(Inew)];
    
    ind = ind + 1;
end

dim_Hilbert = length(state_hash);
if dim_Hilbert > 5e3
    error('dimension = %d,  too large (for N1 = [%d, %d])', dim_Hilbert, N1(1), N1(2));
end

% just FYI for debugging purpose
% fprintf('dim_ResonantSector = %d, dim_Hilbert = %d\n', dim_ResonantSector, dim_Hilbert);

% rearrange state_hash ordering so the resonant sector are at the top
[~, Iresonant, ~] = intersect(state_hash, resonantSector);
Ioffres = setdiff(1:dim_Hilbert, Iresonant);
state_hash = [state_hash(Iresonant); state_hash(Ioffres)];

%% Now fill in the matrix elements of the Hamiltonian

Ham = zeros(dim_Hilbert);

for ind = 1:dim_ResonantSector
    [~,ia,ib] = intersect(state_hash, state_hashes{ind});
    mat_ind = find(state_hash==resonantSector(ind));
    Ham((ia-1)*(dim_Hilbert+1)+1) = dEvecs{ind}(ib); % trick to set the diagonal
    Ham(mat_ind, ia) = Vddvecs{ind}(ib);
    Ham(ia, mat_ind) = Vddvecs{ind}(ib)';
end

%% diagonalize resonant sector to get asymptotic vectors and other info
info.HamResSect = Ham(1:dim_ResonantSector, 1:dim_ResonantSector)*R^3;
[V, ~] = eig(info.HamResSect);
target_vec = sparse(dim_Hilbert, dim_ResonantSector);
target_vec(1:dim_ResonantSector, 1:dim_ResonantSector) = V;

info.overlap_vectors = target_vec;


%% Finally, diagonalize the full Hamiltonian to get energies and overlaps
    
[V, D] = eig(Ham);
overlaps = V'*target_vec;
[~, I] = max(abs(overlaps(:,1))); % find the state closest to the first state in the resonant sector
energies = diag(D);
y = energies(I);

info.state_labels = getStateLabel(state_hash);

% === END OF MAIN FUNCTION ===
%% now for helper functions/subroutines

    function [dE, Vdd, state_hash_ind] = getInteractions(N1, L1, J1, M1)
        %get a list of states that interact with the two-body state labeled
        %by input, outputing dE (Forster defects) and Vdd (dipole interaction
        %matrix element), and state_hash_ind (corresponding state indices)
        
        N2=cell(2,1); L2=N2; J2=N2; M2=N2;
        dE=N2;dp=N2;
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

            if abs(diff(N1)) > N2s_distance
                N2s = union(N2s, N1(3-r)+(-N2s_distance:N2s_distance));
            end

            % The "forbidden" are: (i) dJ=0,+1,-1   (ii) J2 not negative, and can support M2 (iii) L2 not negative
            % JDT: Now remove the states that are not physical (in terms of angular
            % momentum quantum numbers only)
            Iforbid=(abs(J2s-J1(r))>1 | J2s<abs(M2s) | L2s<0);L2s(Iforbid)=[];J2s(Iforbid)=[];M2s(Iforbid)=[];

            % JDT: now construct arrays of angular momentum quantum numbers for
            % each principal quantum number
            [~,L2s]=meshgrid(N2s,L2s); [~,J2s]=meshgrid(N2s,J2s); [N2s,M2s]=meshgrid(N2s,M2s);
            N2s=N2s(:);J2s=J2s(:);L2s=L2s(:);M2s=M2s(:); % JDT: Flatten into one-dim array

            % JDT: now remove unphysical states where L>N
            IlargeL=(L2s>=N2s);N2s(IlargeL)=[];L2s(IlargeL)=[];J2s(IlargeL)=[];M2s(IlargeL)=[];

            N2{r} = N2s; L2{r} = L2s; J2{r} = J2s; M2{r} = M2s;

            % JDT: get angular part of matrix element
            reduced_angular=reduced_LJM(L1(r),J1(r),M1(r),L2s,J2s,M2s);

            % JDT: get energy difference in AU
            dE{r}=(nstar(N1(r),L1(r),J1(r)).^(-2)-nstar(N2s,L2s,J2s).^(-2))/2/ReducedMassFactor; % To convert to a.u. (factor 2 because the Hartree energy is 2*Ry);
            N1s=N2s*0+N1(r);
            if any(N1s > 150 | N2s > 150)
                error('large unsupported N');
            end
            Rdipole12=Rdipole_table(N1s,L1(r),J1(r),N2s,L2s,J2s);

            % JDT: calculate total matrix elements, and save magnetic sublevels
            dp{r}=(reduced_angular.*Rdipole12);
        end
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
        geom_factor_free_space = I1p.*I2m+I1m.*I2p+(I1z.*I2z)*(1-3*cos(theta)^2)-...
            (3/2)*sin(theta)^2*(I1p.*I2p+I1p.*I2m+I1m.*I2p+I1m.*I2m)-...
            (3/sqrt(2))*sin(theta)*cos(theta)*(I1p.*I2z+I1m.*I2z+I1z.*I2p+I1z.*I2m); % ref. [9]
        Vdd=geom_factor_free_space.*dp./R^3;
        
        [N2{1} N2{2}] = meshgrid(N2{1}, N2{2}); 
        [L2{1} L2{2}] = meshgrid(L2{1}, L2{2}); 
        [J2{1} J2{2}] = meshgrid(J2{1}, J2{2}); 
        [M2{1} M2{2}] = meshgrid(M2{1}, M2{2}); 
        state_hash_ind = getHashInd([N2{1}(:), N2{2}(:), L2{1}(:), L2{2}(:), J2{1}(:), J2{2}(:), M2{1}(:), M2{2}(:)]);
        dE = dE(:);
        Vdd = Vdd(:);
    end
    
    function hash = getHashInd(statelabel)
    % statelabel must be in form of [n1, n2, l1, l2, j1, j2, m1, m2]
        statelabel(:, 7:8) = statelabel(:, 7:8) + statelabel(:, 5:6);
        statelabel(:, 5:6) = statelabel(:, 5:6) + 1/2;
        hash = statelabel*hashmultipliers;
    end

    function statelabel = getStateLabel(hash)
    % given hash index (indices), return [n1, n2, l1, l2, j1, j2, m1, m2]
        [meshhash, meshprimes] = meshgrid(hash, hashprimes);
        statelabel = mod(meshhash, meshprimes)';
        statelabel(:, 5:6) = statelabel(:, 5:6) - 1/2;
        statelabel(:, 7:8) = statelabel(:, 7:8) - statelabel(:, 5:6);
    end

end

%% References
% [9] Reinhard, ..., Raithel PRA 75, 032712 (2007).
