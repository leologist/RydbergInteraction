function y=angerj(nu,z,D)
% Anger function of real order nu, and its first derivative
% D=0: J_nu(Z) 
% D=1: J_nu'(Z)
%
% Vectorization wrap.
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.

if isscalar(nu)
    y=angerj_scalar_nu(nu,z,D);
else
    y=z*0;
    for ind=1:length(y(:))
        y(ind)=angerj_scalar_nu(nu(ind),z(ind),D);
    end
end
