%参量不稳定性色散关系矩阵the dispersion matrix of PDI
function y = PIM(w)
w1 = w-1;
w2 = w+1;
global E0x E0y E0z wce w0 n0x n0z;
global nx ny nz nsqu n1x n1y n1z n1squ n2x n2y n2z n2squ;
rn = [nx ny nz];
rn1 = [n1x n1y n1z];
rn2 = [n2x n2y n2z];
% 扰动速度
b = -30./w0.*0.0048./8.19.*[E0x E0y E0z].';
w0x = 30./w0.*0.0048./8.19.*(-n0z.*E0y);
w0y = 30./w0.*0.0048./8.19.*(n0z.*E0x-n0x.*E0z);
w0z = 30./w0.*0.0048./8.19.*(n0x.*E0y);
A = [0 wce+w0z -w0y
    -(wce+w0z) 0 w0x
    w0y -w0x -1i];
v0 = A\b;

% 电磁KEM
d = D(w);%linear dielectric tensor matrix of the low frequency wave
d1 = D1(w);%linear dielectric tensor matrix of the lower sideband wave

s1 = SigNL1(w);%quasi-linear dielectric tensor matrix of the lower sideband wave

s3 = SigNLp1(w);%quasi-linear dielectric tensor matrix of the low frequency wave

sls1 = SigNLls1(w);%nonlinear dielectric tensor matrix of the lower sideband wave
% y=d-s3*inv(d1)*s1-s4*inv(d2)*s2;
y = [d(1,1) d(1,2) d(1,3) -s3(1,1) -s3(1,2) -s3(1,3)
    d(2,1) d(2,2) d(2,3) -s3(2,1) -s3(2,2) -s3(2,3)
    d(3,1) d(3,2) d(3,3) -s3(3,1) -s3(3,2) -s3(3,3)
    s1(1,1) s1(1,2) s1(1,3) sls1(1,1)-d1(1,1) sls1(1,2)-d1(1,2) sls1(1,3)-d1(1,3)
    s1(2,1) s1(2,2) s1(2,3) sls1(2,1)-d1(2,1) sls1(2,2)-d1(2,2) sls1(2,3)-d1(2,3)
    s1(3,1) s1(3,2) s1(3,3) sls1(3,1)-d1(3,1) sls1(3,2)-d1(3,2) sls1(3,3)-d1(3,3)];
%{
% 静电KES
d = (rn./norm(rn))*D(w)*(rn./norm(rn)).';%linear relative dielectric constants of the low frequency wave
d1 = (rn1./norm(rn1))*D1(w)*(rn1./norm(rn1)).';%linear relative dielectric constants of the lower sideband wave

s1 = (rn1./norm(rn1))*SigNL1(w)*(rn./norm(rn1)).';%quasi-linear susceptibility of the lower sideband wave

s3 = (rn./norm(rn))*SigNLp1(w)*(rn1./norm(rn)).';%quasi-linear susceptibility of the low frequency wave

sls1 = (rn1./norm(rn1))*SigNLls1(w)*(rn1./norm(rn1)).';%nonlinear susceptibility of the lower sideband wave

y = [d -s3
    s1 sls1-d1];
%}
end