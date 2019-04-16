function w = w6j(j1,j2,j3,j4,j5,j6)

% W = W6J(j1,j2,j3,j4,j5,j6)
%
% Calculates Wigner 6-j Symbol using Racah W-coefficients
%       { j1 j2 j3 }
%       { j4 j5 j6 }
%
% (j1,j2,j3),(j1,j5,j6),(j4,j2,j6) and (j4,j3,j5) must satisfy triangle
% condition |j2 - j3|<=j1<=j2 + j3
%
% All square roots use logarithms to cope with large factorials
%
% Equations (7.25-7) in 'Angular Momentum:Techniques in Quantum Mechanics'
% by V. Devanathan, Kluwer Academic Publishers (1999).
%
%  Vectorized by: O. Firstenberg, Harvard University HQOC ; MIT, 2012.
%
%  Based on original code by:
%    J. Pritchard Durham University 2009
%    [downloaded from: http://massey.dur.ac.uk/jdp/code.html, Jan. 2012]


if triangle(j1,j2,j3)==1
    error(sprintf('j1,j2,j3 must satisfy triangle relation:\n\n\t\t|j2 - j3|<=j1<=j2 + j3'));
elseif triangle(j1,j5,j6)==1
    error(sprintf('j1,j5,j6 must satisfy triangle relation:\n\n\t\t|j5 - j6|<=j1<=j5 + j6'));
elseif triangle(j4,j2,j6)==1
    error(sprintf('j4,j2,j6 must satisfy triangle relation:\n\n\t\t|j2 - j6|<=j4<=j2 + j6'))
elseif triangle(j4,j3,j5)==1
    error(sprintf('j4,j3,j4 must satisfy triangle relation:\n\n\t\t|j3 - j5|<=j4<=j3 + j5'));
else
    wracah=j1*0;wracah(:)=racah(j1(:),j2(:),j5(:),j4(:),j3(:),j6(:));
    w = (-1).^(j1+j2+j4+j5).*wracah;
end

%Racah Coefficient
function W = racah(a,b,c,d,e,f)
W = tri(a,b,e).*tri(c,d,e).*tri(a,c,f).*tri(b,d,f);

x=(-1:max([(a+b+c+d),(a+d+e+f),(b+c+e+f)]));
[X,a]=meshgrid(x,a);[X,b]=meshgrid(x,b);[X,c]=meshgrid(x,c);
[X,d]=meshgrid(x,d);[X,e]=meshgrid(x,e);[X,f]=meshgrid(x,f);
I= X>=(a+b+e) & X>=(c+d+e) & X>=(a+c+f) & X>=(b+d+f) & X<=(a+b+c+d) & X<=(a+d+e+f) & X<=(b+c+e+f);
xsum=X*0;
x=X(I);a=a(I);b=b(I);c=c(I);d=d(I);e=e(I);f=f(I);
xsum(I)=(-1).^(x+a+b+c+d)...
        .*exp(lgf(x+1)-lgf(x-a-b-e)-lgf(x-c-d-e)...
        -lgf(x-a-c-f)-lgf(x-b-d-f)-lgf(a+b+c+d-x)...
        -lgf(a+d+e+f-x)-lgf(b+c+e+f-x));
xsum=sum(xsum,2);
W = W.*xsum; 

%Triangle Coefficient
function t=tri(a,b,c)
% t = sqrt(f(a+b-c)*f(a-b+c)*f(-a+b+c)/f(a+b+c+1));
t = exp(0.5*(lgf(a+b-c)+lgf(a-b+c)...
    +lgf(-a+b+c)-lgf(a+b+c+1)));

%Triangle relations test |a-b|<=c<=a+b
function test=triangle(a,b,c)
test=any((a<abs(b-c))|(a>(b+c)));