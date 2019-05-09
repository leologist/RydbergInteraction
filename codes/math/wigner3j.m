function w=w3j(j1,m1,j2,m2,j3,m3)

% W = W3J(j1,m1,j2,m2,j3,m3)
%
% j1,j2,j3 must satisfy triangle condition |j2 - j3|<=j1<=j2 + j3
%
% Wigner 3-j symbol is evaluated using equation found in 'Angular 
% Momentum: An Illustrated guide to Rotational Symmetries for
% Physical Systems', W. J. Thompson
%
%  Vectorized by: O. Firstenberg, Harvard University HQOC ; MIT, 2012.
% 
%  Based on original code by:
%    J. Pritchard Durham University 2009
%    [downloaded from: http://massey.dur.ac.uk/jdp/code.html, Jan. 2012]


%Check Triangular relation
if any((j3<abs(j1-j2))|(j3>(j1+j2)))
    error(sprintf('Addition of angular momentum requires triangle relation\n\t|j1-j2|<=j3<=j1+j2'));
%Evaluate w3j
else
     wksum=m1*0;
     wksum(:)=ksum(j1(:),m1(:),j2(:),m2(:),j3(:),m3(:));
     w=(-1).^(j1-j2-m3).*...
         exp(0.5*(lgf(j3+j1-j2)+lgf(j3-j1+j2)+lgf(j1+j2-j3)...
         +lgf(j3-m3)+lgf(j3+m3) - lgf(j1+j2+j3+1)...
         -lgf(j1-m1)-lgf(j1+m1)-lgf(j2-m2)-lgf(j2+m2))).*wksum;
    w(m1+m2+m3~=0)=0;
    w(j2==0)=1;
end

%Summation performed for all values of k which give non-negative fs
function s=ksum(j1,m1,j2,m2,j3,m3)
k=(0:max([max(j3-j1+j2) max(j3-m3)]));
[K,j1]=meshgrid(k,j1);[K,j2]=meshgrid(k,j2);[K,j3]=meshgrid(k,j3);
[K,m1]=meshgrid(k,m1);[K,m2]=meshgrid(k,m2);[K,m3]=meshgrid(k,m3);
I=(K>=-j1+j2-m3) & (K<=j3-j1+j2) & (K<=j3-m3);
s=K*0;
K=K(I);j1=j1(I);j2=j2(I);j3=j3(I);m1=m1(I);m2=m2(I);m3=m3(I);
s(I)=(-1).^(K+j2+m2).*exp(lgf(j2+j3+m1-K)+lgf(j1-m1+K)-lgf(K)-lgf(j3-j1+j2-K)-lgf(j3-m3-K)-lgf(K+j1-j2+m3));
s=sum(s,2);