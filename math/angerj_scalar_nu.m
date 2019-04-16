function y=angerj_scalar_nu(nu,z,D)
% Anger function of real order nu, and its first derivative
%
% D=0: J_nu(Z) 
% D=1: J_nu'(Z)
%
% O. Firstenberg, Harvard University HQOC ; MIT, 2012.

if round(nu)==nu
    switch D
        case 0, y=besselj(nu,z);
        case 1, y=(besselj(nu-1,z)-besselj(nu+1,z))/2;
        otherwise, error('Missing implementation for this derivative order');
    end
else
    switch D
        case 0
            y=   pfq(1,[  1 - nu/2, nu/2 + 1  ], -(z.^2/4),0)./(       nu        )-...
                z.*pfq(1,[3/2 - nu/2, nu/2 + 3/2], -(z.^2/4),0)./((nu- 1).*(nu + 1));
        case 1
            y=   pfq(1,[  1 - nu/2, nu/2 + 1  ], -(z.^2/4),1)./(       nu        ).*(-z/2)-...
                pfq(1,[3/2 - nu/2, nu/2 + 3/2], -(z.^2/4),0)./((nu- 1).*(nu + 1))-...
                z.*pfq(1,[3/2 - nu/2, nu/2 + 3/2], -(z.^2/4),1)./((nu- 1).*(nu + 1)).*(-z/2);
        otherwise
            error('Missing implementation for this derivative order');
    end
    y=y.*2.*sin((pi*nu)/2).*cos((pi*nu)/2)/pi;
end
