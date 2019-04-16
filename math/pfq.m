function [y, tt, nterms] = pfq(a,b,z,d)
%PFQ Generalized Hypergeometric Function pFq
%   [Y,TT,NTERMS] = PFQ(A,B,Z,D) returns the Generalized Hypergeometric
%   Function pFq with numerator parameters A, denominator parameters B,
%   evaluated at the values in Z, with respect to derivatives D.
%
%   PFQ(A,B,Z) returns the values evaluated at Z
%
%   Inputs:
%   A : real or complex valued vector (or empty [])
%   B : real or complex valued vector (or empty [])
%   Z : real or complex valued M-by-N-by-... array
%   D : scalar or vector of real non-negative integers
%
%   Outputs:
%   Y      : Evaluation of pFq at values in Z, with respect to derivative D
%   TT     : Estimated Truth Table on values given
%              +1 : Terms were convergent, answer valid close to precision
%               0 : Uncharacterized, use answers with caution
%              -1 : Terms were divergent, use answers with severe caution
%   NTERMS : Number of terms used in the calculation for each output
%
%   If Z is M-by-N-by-... array and D is a vector 1-by-Q then the outputs
%   will be an array of size M-by-N-by-...-by-Q
%
%   Examples:
%   %0F1: Confluent Hypergeometric Function
%   pfq( [] , 1 , 1.5 )
%   pfq( [] , 1+1i , 1.5-1i )
%   pfq( [] , 1 , -2 )
%   %1F1: Kummer Confluent Hypergeometric Function
%   pfq( 1 , 2 , 3 )
%   pfq( 2+1i , 2, 0.5 )
%   pfq( 10 , 1/3 , -1 )
%   %2F1: Hypergeometric Function
%   pfq( [1/3 1/3] , 2/3 , 1/2 )
%   pfq( [2+1i -1i] , 3/4 , 0.5-0.5i )
%   %PFQ: Generalized Hypergeometric Function
%   pfq( [1 1] , [3 3 3] , 2 )
%   pfq( [1i 1i 1i] , [2 2 2] , -1i )
%   pfq( [1 2 3 4] , [5 6 7] , [.1 .3 .5] )
%
%   Vectorization with derivatives:
%   Evaluate pFq above, at [.1 .3 .5], but with respect to the 0th, 1st and
%   2nd derivative. Show all outputs
%   [y , tt , nterms] = pfq( [1 2 3 4] , [5 6 7] , [.1 .3 .5] , 0:2 )
%
%   Finite degree polynomial, with derivatives
%   x = [.1 .3 .5];
%   y = pfq( [1 -2] , -2 , x , [0:4] )
%   [1+x+x.^2; 1+2.*x] %Theoretical answer
%

%
%   Mike Sheppard
%   MIT Lincoln Laboratory
%   michael.sheppard@ll.mit.edu
%   Last Modified 25-Jul-2011

%
%   Comments, suggestions, or improvements are welcome
%


%Check inputs
if nargin < 3
    error('pfq:TooFewInputs',...
        'Requires at least three input argument.');
end
if nargin == 3
    d=0; %Default zero'th derivative
end
if ~isvector(a) && ~isempty(a)
    error('pfq:A must be either a vector or empty');
end
if ~isvector(b) && ~isempty(b)
    error('pfq:B must be either a vector or empty');
end
if ~(isvector(d) && all(d==round(d)) && all(d>=0))
    error('pfq:D must be either a scalar or vector of real non-negative integers');
end


%Tolerance
tol_exp=50;
%Tolerance Exponent, base 10
%   Terms are added until the norm of each term, for all inputs and
%   derivatives, are ...
%   1. Less than 10^(-TOL_EXP) in absolute value; converging
%   2. Greater than 10^(TOL_EXP) in absolute value; diverging
%   Note: It may be oscillatory and an exact answer exists, however this
%   will stop when the terms exceed 10^(TOL_EXP) in absolute value and
%   yields -1 for Estimated Truth Table, above.

%Maximum number of steps
maxsteps=1e4;
    



%
%   The definition of pFq by a formal power series is:
%   pFq = sum( [C(a,n) / C(b,n)] * [z^n/n!] , n=0, ... , infinity )
%      with
%   C(v,n)=prod( [gamma(v(j)+n) / gamma(v(j)] , j=1, ... ,length(v) )
%
%   The derivatives are calculated equivalently, with [C(a,n)/C(b,n)]
%   term remaining the same and using d'th derivative of
%   z^k/k! -> z^(k-d) / (k-d)!
%
%   Terms are added until all the terms, for all inputs and
%   derivatives, exceed TOL given below; or if expansion terminates
%   with a finite polynomial or maximum allowed number of steps
%


% INITIALIZE PARAMETERS
%Make A and B row vectors for consistency
a=(a(:)).'; b=(b(:)).'; %non-conjugate transpose of column vector
%Keep track of which A's and B's are negative integers
anegint=(a<0)&(a==round(a)); %Negative integer values of A
bnegint=(b<0)&(b==round(b)); %Negative integer values of B
%If any A's is a negative integer a polynomial is produce,
%find stopping point
if any(anegint)
    minnega=min(-a(anegint)); %stopping point
else
    minnega=Inf; %no stopping
end

sz=size(z); ld=length(d); z=z(:); lnz=log(z); lz=length(lnz);
y=zeros([lz ld]); tol_high=10^abs(tol_exp); tol_low=tol_high^-1;
n=0; keepgoing=1;
indxz=true(lz,ld); %Logical index of which Z, and which derivatives, to still work on
tt=zeros(lz,ld); %Estimated Truth Table
nterms=ones(lz,ld); %Number of terms, minimum number is one


while keepgoing
    
    %Compute the log of product of terms
    %C(v,n)=prod( [gamma(v(j)+n) / gamma(v(j)] , j=1, ... ,length(v) )
    %Returns log of product, by summing the individual log of terms
    %Subfunction watches out for negative integers and does those
    %separately
    lnnum=lnprodterms(a,anegint,n);
    lnden=lnprodterms(b,bnegint,n);
    term=zeros([lz ld]); %set to zero
    
    %Loop through each derivative given
    for dn=1:ld
        %Use only Z values that need more computation, for this derivative
        indxzd=indxz(:,dn); z_temp=z(indxzd); findxzd=find(indxzd);
        if ~isempty(findxzd)
            %For each derivative, check to see if term in original expansion is
            %valid part of derivative. The indexing is for the actual function;
            %and does not change for derivatives, but coefficient and exponent
            %do change accordingly
            %d'th derivative of (coef) * z^n / n! -> (coef) * z^(n-d) / (n-d)!
            if n>=d(dn)
                %Rewrite in log space, using d'th derivative
                lnzn=(n-d(dn)).*lnz(indxzd); lnfact=lngamma(n-d(dn)+1);
                %Natural log of term for all z for specific derivative
                lnterm=lnnum-lnden-lnfact+lnzn;
                %Check special case if z==0. If so, and if n==d(dn) then term
                %is coef as computed, else zero. If z~=0 then coef is correct.
                ztz=(z_temp==0); indx0=findxzd(ztz); indx1=findxzd(~ztz);
                %keyboard
                term(indx0,dn)=(n==d(dn)).*exp(lnnum-lnden-lnfact);
                term(indx1,dn)=exp(lnterm(~ztz));
                %Record number of steps
                %Term of original series - #derivative + one
                %as started with zero'th term (z^0)
                nterms(findxzd,dn)=n-d(dn)+1;
                %Terminate if all values are outside range of tolerance;
                %either diverging or converging, keep all those within range to
                %keep on going with terms
                normM=abs(term(indxzd,dn));
                lo_tol=(normM<tol_low); hi_tol=(normM(:,1)>tol_high);
                if any(lo_tol | hi_tol)
                    %Add either +1 or -1 to Estimated Truth Table
                    tt(findxzd(lo_tol),dn)=1; tt(findxzd(hi_tol),dn)=-1;
                    %Delete those points from further calculations
                    indxz(findxzd(lo_tol | hi_tol),dn)=false;
                end
            else
                %If term in original expansion is not a valid part of the
                %derivative then just add zero
                term(:,dn)=0;
            end
        end  %~isempty(findxzd)
    end
    
    %Add the new term
    y=y+term;
    
    %Break out of loop if n is equal to minimum of any A's that are
    %negative integers (finite degree polynomial)
    if n==minnega
        keepgoing=0;
    end
    
    %If all new terms for all values exceed the tolerance; stop.
    if all(indxz(:)==false)
        keepgoing=0; %Stop, all values exceed tolerance
    else
        n=n+1; %Keep going to the next term
    end
    
    %Break out of loop after MAXSTEPS terms, if haven't already
    %This is usually the case when z is close to 1
    if n>maxsteps
        keepgoing=0;
        %Record all remaining ones as uncharacterized, 
        %and maxsteps number of steps. 
        nterms(indxz)=maxsteps;
        tt(indxz)=0; %uncharacterized
    end
    
    
end


%For values of z==0, results are accurate for any derivative
if any(z==0)
    k0=(z==0); tt(k0,:)=1; 
end

%For finite polynomials the results are accurate for any derivative
%except if output exceeded given tolerance
if isfinite(minnega)
    tt(tt~=-1)=1; %keep -1 the same, change all others to 1
end

%For those with MAXSTEPS, change TT to 0 for uncharacterized
tt(nterms==maxsteps)=0;

%Reshape final result to proper dimensions from input
y=reshape(y(:),[sz ld]);
tt=reshape(tt(:),[sz ld]);
nterms=reshape(nterms(:),[sz ld]);
if isvector(z)
    if sz(1)==1 && ld>1
        y=(squeeze(y))';
        tt=(squeeze(tt))';
        nterms=(squeeze(nterms))';
    else
        y=squeeze(y);
        tt=squeeze(tt);
        nterms=squeeze(nterms);
    end
end


%Fix round-off
lim=100*eps;
%Make real part zero, keep imaginary
kr=abs(real(y))<lim; if any(kr), y(kr)=imag(y(kr)); end
%Make imag part zero, keep real part
ki=abs(imag(y))<lim; if any(ki), y(ki)=real(y(ki)); end


return




function lnprod=lnprodterms(v,vnegint,n)
%Given vector A, and term n compute
%Log[ (v_1)(n) ... (v_t)(n) ] where
%(v)(n)= (v)*(v+1)*(v+2)*...*(v+n-1) [Pochhammer symbol]
%Where vnegint is logical index of negative integers within V

%Step 1: Compute non-negative integers (valid input for lngamma)
lnprod=sum(lngamma(v(~vnegint)+n)-lngamma(v(~vnegint)));

%Step 2: Add the remaining terms that use negative integers
if (any(vnegint))&&(n>0)  %Exclude n==0 term
    A=v(vnegint); A=A(:).'; B=0:n-1; %Correct orientation for bsxfun
    sumv=bsxfun(@plus,A,B);
    %Trivial Example: A=[4 7 16].', B=0:3, sumv=[A A A A]+[B;B;B]
    lnprod=lnprod+sum(log(sumv(:))); %sum of logs -> log of products
end

return






%THE CODE BELOW IS FROM
%http://www.mathworks.com/matlabcentral/fileexchange/ ...
%  ... 978-special-functions-math-library/content/gammaln.m
%
%Modifications:
%   1. Renamed to lngamma so as not to override built-in gamma or gammaln
%   2. Replaced i with 1i
% -Mike Sheppard (7/20/2011)
%

function [f] = lngamma(z)
% GAMMALOG  Natural Log of the Gamma function valid in the entire complex plane.
%           This routine uses an excellent Lanczos series approximation
%           for the complex ln(Gamma) function.
%
%usage: [f] = gammaln(z)
%             z may be complex and of any size.
%             Also  n! = prod(1:n) = exp(gammalog(n+1))
%
%tested under version 5.3.1
%
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996
%
%see also:   GAMMA GAMMALN GAMMAINC PSI
%see also:   mhelp GAMMA
%see also:   mhelp lnGAMMA

%Paul Godfrey
%pgodfrey@conexant.com
%07-13-01


siz = size(z);
z=z(:);
zz=z;

f = 0.*z; % reserve space in advance

p=find(real(z)<0);
if ~isempty(p)
    z(p)=-z(p);
end

%Lanczos approximation for the complex plane

g=607/128; % best results when 4<=g<=5

c = [  0.99999999999999709182;
    57.156235665862923517;
    -59.597960355475491248;
    14.136097974741747174;
    -0.49191381609762019978;
    .33994649984811888699e-4;
    .46523628927048575665e-4;
    -.98374475304879564677e-4;
    .15808870322491248884e-3;
    -.21026444172410488319e-3;
    .21743961811521264320e-3;
    -.16431810653676389022e-3;
    .84418223983852743293e-4;
    -.26190838401581408670e-4;
    .36899182659531622704e-5];

s=0;
for k=size(c,1):-1:2
    s=s+c(k)./(z+(k-2));
end

zg=z+g-0.5;
s2pi= 0.9189385332046727417803297;

f=(s2pi + log(c(1)+s)) - zg + (z-0.5).*log(zg);

f(z==1 | z==2) = 0.0;

if ~isempty(p)
    lpi= 1.14472988584940017414342735 + 1i*pi;
    f(p)=lpi-log(zz(p))-f(p)-log(sin(pi*zz(p)));
end

p=find(round(zz)==zz & imag(zz)==0 & real(zz)<=0);
if ~isempty(p)
    f(p)=Inf;
end

f=reshape(f,siz);

return
