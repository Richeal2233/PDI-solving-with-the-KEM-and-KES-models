%主程序
clear all;
clc;
fidin=fopen('indata.txt');%JET
fidout=fopen('inmatlab.txt','w');%initialized parameters
while ~feof(fidin)
    tline=fgetl(fidin);
    if double(tline(1))>=48&&double(tline(1))<=57
        fprintf(fidout,'%s\n\n',tline);
        continue
    end
end
fclose(fidout);
out=importdata('inmatlab.txt');
global w0 n0z de1 n1z;
Bs=out(1);%static magnetic field B_T in T
P0=out(2);%pump power P_0 in MW
Te=out(3);%electron temperature T_e in eV
Ti=out(4);%ion temperature in T_i eV
ne=out(5);%electron density in n_e 1e18 m^-3; equal to ion density
w0=out(6);%pump frequency f_0 in GHz(without 2*pi)
n0z=out(7);%parallel refractive index of the pump wave n_{0z}
de1=out(8)./180.*pi;%angle between the perpendicular refractive index of the pump wave and the lower sideband wave in circular measure, default pi/2(90°)
n1z=out(9);%parallel refractive index of the lower sideband wave n_{1z}

nperstart=7.14;%perpendicular refractive index of the low frequency wave n_{LF⊥}
nperend=7.17;
nper_N = 10;
nperdata=linspace(nperstart,nperend,nper_N);

%等离子体参数的计算calculation of the plasma parameters
global wpe wpi wce wci;
global n0squ n0x nz n2z;
nz=n0z+n1z;
n2z=n0z+nz;
mratio=1836;
% mratio = 2.*1836; % D plasma
w0=w0.*2.*pi;
wpe = 56.41.*sqrt(ne)./w0;%electron plasma frequency
% wpe = sqrt(ne).*4.8.*sqrt(40.*pi./0.91)./w0;
wpi=wpe./sqrt(mratio);%ion plasma frequency
wce = 1600.*Bs./9.11./w0;%electron cyclotron frequency
% wce = Bs.*160./0.91./w0;
wci=wce./mratio;%ion cyclotron frequency
%w0=w0./w0;

% 电磁
eps0xx=1+wpe.^2./wce.^2-wpi.^2;
eps0xy=-1i.*wpe.*wpe./wce;
eps0zz=1-wpe.^2;
A0=eps0xx;
B0=(eps0xx+eps0zz).*(n0z.^2-eps0xx)-eps0xy.^2;
C0=eps0zz.*(eps0xy.^2+(n0z.^2-eps0xx).^2);
n0x=sqrt((-B0+sqrt(B0.^2-4.*A0.*C0))./2./A0);
n0squ=n0x.^2+n0z.^2;

% 静电
%{
n0squ=n0z.^2./((1+wpe.^2./wce.^2)./wpi.^2-1).*mratio;
n0x=sqrt(n0squ-n0z.^2);
%}

%热速度
global cs ve vi;
cs = sqrt(Te.*1.6./1.67)./3.*1e-4;
% cs = sqrt(Te.*1.6./1.67./2)./3.*1e-4; % D plasma
vi = sqrt(2.*Ti./Te).*cs;
ve = sqrt(2.*mratio).*cs;
% ve = sqrt(2.*Te.*1.6./9.1)./300;
% vi = ve.*sqrt(Ti./Te./mratio);

%泵浦波电场的计算
global E0x E0y E0z u;
% JET
LyLz = 0.348*0.884;
% EAST
% LyLz = 0.472*0.546;
% C-Mod
% LyLz = 0.24.*0.1525;
% FTU
% LyLz = 1;
e0 = 8.8419e-12;
% e0 = 8.8542e-12;
ratio1=2.*(1+wpe.^2./wce.^2).*e0;
ratio2=3e+8./(1+mratio.*n0z.^2./n0squ).*mratio.*n0z.^2.*n0x./n0squ.^2;
E0squ=8.*pi.*P0.*1e+6./LyLz./ratio1./ratio2;

% 电磁
E0x=sqrt(E0squ./(1+(eps0xy./(eps0xx-n0squ)).^2+(n0x.*n0z./(eps0zz-n0x.^2)).^2));
E0y=E0x.*eps0xy./(eps0xx-n0squ);
E0z=-E0x.*n0x.*n0z./(eps0zz-n0x.^2);

% 静电
%{
E0x=sqrt(E0squ./(1+n0z.^2./n0x.^2));
E0y=0;
E0z=sqrt(E0squ-E0x.^2);
%}

u = E0x./Bs./3e+08;

% 电场强度由国际单位制换算为高斯单位制
E0x = 1./3e+04.*E0x;
E0y = 1./3e+04.*E0y;
E0z = 1./3e+04.*E0z;

%解参量不稳定性色散关系的主程序
lb=[7.375,0.087];%lower bound
ub=[13.375,2.087];%upper bound
PSoptions = optimoptions('patternsearch', 'TolMesh', 1e-3, 'InitialMeshSize', 1, 'Display', 'off');%options for patternsearch
wiv=[10.3750,1.0870];%normlized by 1e6, easy to operate

global nper n1per n1x n1y n1squ nx de ny nsqu n2x n2y n2per de2 n2squ;


% 解第一个根
nper=nperdata(1);
n1per=-n0x.*cos(de1)+sqrt((n0x.*cos(de1)).^2-(n0x.^2-nper.^2));
n1x=n1per.*cos(de1);
n1y=n1per.*sin(de1);
n1squ=n1per.^2+n1z.^2;
nx=n0x+n1x;
ny=n1y;
de=asin(ny./nper);
nsqu=nx.^2+ny.^2+nz.^2;
n2x=n0x+nx;
n2y=ny;
n2per=sqrt(n2x.^2+n2y.^2);
de2=asin(n2y./n2per);
n2squ=n2x.^2+n2y.^2+n2z.^2;

[wans,fval,exitflag] = patternsearch(@(x) fw(x),wiv,[],[],[],[],lb,ub,[],PSoptions);%solving process

fprintf('exitflag=%f',exitflag);
fprintf('fval=%f',fval);
if(exitflag<0)
    wans(1) = 0;
end

wansre = wans(1)*1e-3;%(turn to normlized by 1e9, keep consistent)
wansim = wans(2)*1e-3;%(turn to normlized by 1e9, keep consistent)
fprintf('w_re=%f',wansre);%real frequency of low frequency wave
fprintf('gamma=%f',wansim);%growth rate of low frequency wave
%{
% 解第二个根,根据第一个根
nper=nperdata(2);
n1per=-n0x.*cos(de1)+sqrt((n0x.*cos(de1)).^2-(n0x.^2-nper.^2));
n1x=n1per.*cos(de1);
n1y=n1per.*sin(de1);
n1squ=n1per.^2+n1z.^2;
nx=n0x+n1x;
ny=n1y;
de=asin(ny./nper);
nsqu=nx.^2+ny.^2+nz.^2;
n2x=n0x+nx;
n2y=ny;
n2per=sqrt(n2x.^2+n2y.^2);
de2=asin(n2y./n2per);
n2squ=n2x.^2+n2y.^2+n2z.^2;
[wans(2),fval,exitflag] = fsolve(@fw,wans(1).*1e+3,options);
% [wans(2),fval,exitflag] = fsolve(@fw,wans(1).*1e+3);
fprintf('exitflag=%f',exitflag);
fprintf('fval=%f',fval);
if(exitflag<0)
    wans(2) = 0;
end
wans(2) = wans(2).*1e-3;
wansre(2) = real(wans(2));
wansim(2) = imag(wans(2));
% 解后面的根
for ii=3:nper_N
    nper=nperdata(ii);
    n1per=-n0x.*cos(de1)+sqrt((n0x.*cos(de1)).^2-(n0x.^2-nper.^2));
    n1x=n1per.*cos(de1);
    n1y=n1per.*sin(de1);
    n1squ=n1per.^2+n1z.^2;
    nx=n0x+n1x;
    ny=n1y;
    de=asin(ny./nper);
    nsqu=nx.^2+ny.^2+nz.^2;
    n2x=n0x+nx;
    n2y=ny;
    n2per=sqrt(n2x.^2+n2y.^2);
    de2=asin(n2y./n2per);
    n2squ=n2x.^2+n2y.^2+n2z.^2;
    if(abs(wans(ii-1))<1e-19)
        [wans(ii),fval,exitflag] = fsolve(@fw,wiv,options);
        %[wans(ii),fval,exitflag] = fsolve(@fw,wiv);
    else
        [wans(ii),fval,exitflag] = fsolve(@fw,(2.*wans(ii-1)-wans(ii-2)).*1e+3,options);
        %[wans(ii),fval,exitflag] = fsolve(@fw,(2.*wans(ii-1)-wans(ii-2)).*1e+3);
    end
    fprintf('exitflag=%f',exitflag);
    fprintf('fval=%f',fval);
    if(exitflag<0)
        wans(ii) = 0;
    end
    wans(ii) = wans(ii).*1e-3;
    wansre(ii) = real(wans(ii));
    wansim(ii) = imag(wans(ii));
end
%}

%{
for ii=1:nper_N
    nper=nperdata(ii);
    n1per=-n0x.*cos(de1)+sqrt((n0x.*cos(de1)).^2-(n0x.^2-nper.^2));
    n1x=n1per.*cos(de1);
    n1y=n1per.*sin(de1);
    n1squ=n1per.^2+n1z.^2;
    nx=n0x+n1x;
    ny=n1y;
    de=asin(ny./nper);
    nsqu=nx.^2+ny.^2+nz.^2;
    n2x=n0x+nx;
    n2y=ny;
    n2per=sqrt(n2x.^2+n2y.^2);
    de2=asin(n2y./n2per);
    n2squ=n2x.^2+n2y.^2+n2z.^2;
    [wans(ii),fval,exitflag] = fsolve(@fw,wiv,options);
    fprintf('exitflag=%f',exitflag);
    fprintf('fval=%f',fval);
    if(exitflag<0)
        wans(ii) = 0;
    end
    wans(ii) = wans(ii).*1e-3;
    wansre(ii) = real(wans(ii));
    wansim(ii) = imag(wans(ii));
end
%}