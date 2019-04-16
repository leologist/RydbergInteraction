function y=scalar_pol(n,L,J,M)
%
% Calculates scalar-polarizability for the state |nLJM>
%
% Based on Refs. [8+10]
%
% Vectorized, with array shape determined by the array n.
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.

global atom_name
switch atom_name
    case 'Rb'
        %S J  M    beta1   beta2
        stark_table=[...
            0 1/2 1/2  2.188   5.486;...
            1 1/2 1/2  2.039  51.456;...
            1 3/2 1/2  2.449  62.011;...
            1 3/2 3/2  1.611  52.948;...
            2 3/2 1/2  2.694  -6.159;...
            2 3/2 3/2  1.725  22.259;...
            2 5/2 1/2  2.770 -12.223;...
            2 5/2 3/2  2.352   1.772;...
            2 5/2 5/2  1.513  29.763;...
            3 5/2 1/2 -1.655   1.612e3;...
            3 5/2 3/2 -1.308   1.350e3;...
            3 5/2 5/2 -0.634   0.826e3;...
            3 7/2 1/2 -1.624   1.623e3;...
            3 7/2 3/2 -1.457   1.478e3;...
            3 7/2 5/2 -1.077   1.188e3;...
            3 7/2 7/2 -0.530   0.753e3];
        stark_table(:,4)=stark_table(:,4)*10^-9;
        stark_table(:,5)=stark_table(:,5)*10^-11;
    otherwise
        error('No scalar-polarizability data for this atom');
end
if isscalar(L), L=n*0+L; end
if isscalar(J), J=n*0+J; end
if isscalar(M), M=n*0+M; end
    
table_row=n*0;
for ind=1:size(stark_table,1)
    table_row((L==stark_table(ind,1) & J==stark_table(ind,2) & M==stark_table(ind,3)))=ind;
end
nstar1=nstar(n,L,J);
beta1=stark_table(table_row,4);
beta2=stark_table(table_row,5);

if sum(abs(size(beta1)-size(n)))>0, beta1=beta1.';beta2=beta2.'; end

y=beta1.*nstar1.^6+beta2.*nstar1.^7;

%% References
% [8] O'Sullivan & Stoiche?, PRA 31, 2718 (1985).
% [10] J. Pritchard, PhD thesis, Durham University (2009).