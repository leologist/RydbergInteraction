function nstar=nstar(n,L,J)
%
% Calculates effective principle quantum numbers n*
%
% Vectorized, with array shape determined by the array n.
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.

global atom_name
switch atom_name
    case 'Rb'
        %state:  S_1/2       P_1/2    P_3/2     D_3/2       D_5/2
        delta0=[3.1311804 2.6548849 2.6416737  1.34809171  1.34646572]; %[ref 2]
        delta2=[0.1784    0.2900    0.2950    -0.60286    -0.59600   ]; %[ref 2]
        
        %state:           F_5/2      F_7/2  L>3
        delta0=[delta0 0.0165192 0.0165437   0   ]; %[ref 3]
        delta2=[delta2    -0.085    -0.086   0   ]; %[ref 3]
    otherwise
        error('No n* data for this atom');
end

if isscalar(L), L=n*0+L; end
if isscalar(J), J=n*0+J; end

table_col=n*0;
table_col((L==0 & J==1/2))=1;
table_col((L==1 & J==1/2))=2;
table_col((L==1 & J==3/2))=3;
table_col((L==2 & J==3/2))=4;
table_col((L==2 & J==5/2))=5;
table_col((L==3 & J==5/2))=6;
table_col((L==3 & J==7/2))=7;
table_col((L>=4 & J==(L+1/2)))=8;
table_col((L>=4 & J==(L-1/2)))=8;

d0=delta0(table_col);
d2=delta2(table_col);
if sum(abs(size(d0)-size(n)))>0, d0=d0.';d2=d2.';end

delta=d0+d2./(n-d0).^2;

nstar=n-delta; 

% References
% [2] Li, Mourachko, Noel, Gallagher, PRA 67, 052502 (2003).
% [3] Han, ..., Gallagher, PRA 74, 054502 (2006).