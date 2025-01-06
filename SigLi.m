%离子 低频 线性
%完成于2015/9/16
function y = SigLi(w)
global wpi wci nz vi de nper;
epsilon = 1e-09;
SigLxx=0;
SigLxy=0;
SigLxz=0;
SigLyx=0;
SigLyy=0;
SigLyz=0;
SigLzx=0;
SigLzy=0;
SigLzz=0;
for n=-0:0
    zn=(w-n.*wci)./nz./vi;
    b=(nper.*vi./wci).^2./2;
    vipar = wci./nper./vi;
    viper=nper.*vi./2./wci;
    bd=(n./b.*besseli(n,b)+besseli(n+1,b)).*exp(-b);
    bd1=((n./b-1).*besseli(n,b)+besseli(n+1,b)).*exp(-b);
    bd2=((2+n.*(n-2.*b-1)./b.^2).*besseli(n,b)-(2+1./b).*besseli(n+1,b)).*exp(-b);
    %xx
    c1=n.^2./b.*besseli(n,b).*exp(-b);
    c2=b.*bd2+bd;
    SigLxx=SigLxx+Z(zn).*(c1.*cos(de).^2+c2.*sin(de).^2);
    %xy
    c3=n.^2./b.*besseli(n,b).*exp(-b)-b.*bd2-bd;
    SigLxy=SigLxy-1i.*Z(zn).*(-n.*bd1+1i.*sin(de).*cos(de).*c3);
    %xz
    c4=1i.*n.*vipar.*besseli(n,b).*exp(-b);
    c5=viper.*bd1;
    SigLxz=SigLxz-2i.*(1+zn.*Z(zn)).*(c4.*cos(de)-c5.*sin(de));
    %yx
    SigLyx=SigLyx+1i.*Z(zn).*(-n.*bd1-1i.*sin(de).*cos(de).*c3);
    %yy
    SigLyy=SigLyy+Z(zn).*(c1.*sin(de).^2+c2.*cos(de).^2);
    %yz
    SigLyz=SigLyz-2i.*(1+zn.*Z(zn)).*(c4.*sin(de)+c5.*cos(de));
    %zx
    SigLzx=SigLzx-2i.*(1+zn.*Z(zn)).*(c4.*cos(de)+c5.*sin(de));
    %zy
    SigLzy=SigLzy-2i.*(1+zn.*Z(zn)).*(c4.*sin(de)-c5.*cos(de));
    %zz
    SigLzz=SigLzz+2.*zn.*besseli(n,b).*exp(-b).*(1+zn.*Z(zn));
end
y = wpi.^2./w./nz./vi.*[SigLxx SigLxy SigLxz
    SigLyx SigLyy SigLyz
    SigLzx SigLzy SigLzz];
yold = y-100;
num = 1;
while max(max(abs(y-yold))) >= epsilon
    yold = y;
    SigLxx=0;
    SigLxy=0;
    SigLxz=0;
    SigLyx=0;
    SigLyy=0;
    SigLyz=0;
    SigLzx=0;
    SigLzy=0;
    SigLzz=0;
    for n=-num:num
        zn=(w-n.*wci)./nz./vi;
        b=(nper.*vi./wci).^2./2;
        vipar = wci./nper./vi;
        viper=nper.*vi./2./wci;
        bd=(n./b.*besseli(n,b)+besseli(n+1,b)).*exp(-b);
        bd1=((n./b-1).*besseli(n,b)+besseli(n+1,b)).*exp(-b);
        bd2=((2+n.*(n-2.*b-1)./b.^2).*besseli(n,b)-(2+1./b).*besseli(n+1,b)).*exp(-b);
        %xx
        c1=n.^2./b.*besseli(n,b).*exp(-b);
        c2=b.*bd2+bd;
        SigLxx=SigLxx+Z(zn).*(c1.*cos(de).^2+c2.*sin(de).^2);
        %xy
        c3=n.^2./b.*besseli(n,b).*exp(-b)-b.*bd2-bd;
        SigLxy=SigLxy-1i.*Z(zn).*(-n.*bd1+1i.*sin(de).*cos(de).*c3);
        %xz
        c4=1i.*n.*vipar.*besseli(n,b).*exp(-b);
        c5=viper.*bd1;
        SigLxz=SigLxz-2i.*(1+zn.*Z(zn)).*(c4.*cos(de)-c5.*sin(de));
        %yx
        SigLyx=SigLyx+1i.*Z(zn).*(-n.*bd1-1i.*sin(de).*cos(de).*c3);
        %yy
        SigLyy=SigLyy+Z(zn).*(c1.*sin(de).^2+c2.*cos(de).^2);
        %yz
        SigLyz=SigLyz-2i.*(1+zn.*Z(zn)).*(c4.*sin(de)+c5.*cos(de));
        %zx
        SigLzx=SigLzx-2i.*(1+zn.*Z(zn)).*(c4.*cos(de)+c5.*sin(de));
        %zy
        SigLzy=SigLzy-2i.*(1+zn.*Z(zn)).*(c4.*sin(de)-c5.*cos(de));
        %zz
        SigLzz=SigLzz+2.*zn.*besseli(n,b).*exp(-b).*(1+zn.*Z(zn));
    end
    y = wpi.^2./w./nz./vi.*[SigLxx SigLxy SigLxz
        SigLyx SigLyy SigLyz
        SigLzx SigLzy SigLzz];
    num = num+1;
end
y = y;