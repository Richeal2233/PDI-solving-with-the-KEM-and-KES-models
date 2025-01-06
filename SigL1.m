%电子 下边频 线性
%完成于2015/9/14
function y = SigL1(w)
w1=w-1;
global wpe wce n1z ve de1 n1per;
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
    zn1=(w1-n.*wce)./n1z./ve;
    b1=(n1per.*ve./wce).^2./2;
    vepar = wce./n1per./ve;
    veper=n1per.*ve./2./wce;
    bd=(n./b1.*besseli(n,b1)+besseli(n+1,b1)).*exp(-b1);
    bd1=((n./b1-1).*besseli(n,b1)+besseli(n+1,b1)).*exp(-b1);
    bd2=((2+n.*(n-2.*b1-1)./b1.^2).*besseli(n,b1)-(2+1./b1).*besseli(n+1,b1)).*exp(-b1);
    %xx
    c1=n.^2./b1.*besseli(n,b1).*exp(-b1);
    c2=b1.*bd2+bd;
    SigLxx=SigLxx+Z(zn1).*(c1.*cos(de1).^2+c2.*sin(de1).^2);
    %xy
    c3=n.^2./b1.*besseli(n,b1).*exp(-b1)-b1.*bd2-bd;
    SigLxy=SigLxy-1i.*Z(zn1).*(n.*bd1+1i.*sin(de1).*cos(de1).*c3);
    %xz
    c4=1i.*n.*vepar.*besseli(n,b1).*exp(-b1);
    c5=veper.*bd1;
    SigLxz=SigLxz-2i.*(1+zn1.*Z(zn1)).*(c4.*cos(de1)+c5.*sin(de1));
    %yx
    SigLyx=SigLyx+1i.*Z(zn1).*(n.*bd1-1i.*sin(de1).*cos(de1).*c3);
    %yy
    SigLyy=SigLyy+Z(zn1).*(c1.*sin(de1).^2+c2.*cos(de1).^2);
    %yz
    SigLyz=SigLyz-2i.*(1+zn1.*Z(zn1)).*(c4.*sin(de1)-c5.*cos(de1));
    %zx
    SigLzx=SigLzx-2i.*(1+zn1.*Z(zn1)).*(c4.*cos(de1)-c5.*sin(de1));
    %zy
    SigLzy=SigLzy-2i.*(1+zn1.*Z(zn1)).*(c4.*sin(de1)+c5.*cos(de1));
    %zz
    SigLzz=SigLzz+2.*zn1.*besseli(n,b1).*exp(-b1).*(1+zn1.*Z(zn1));
end
y = wpe.^2./w1./n1z./ve.*[SigLxx SigLxy SigLxz
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
        zn1=(w1-n.*wce)./n1z./ve;
        b1=(n1per.*ve./wce).^2./2;
        vepar = wce./n1per./ve;
        veper=n1per.*ve./2./wce;
        bd=(n./b1.*besseli(n,b1)+besseli(n+1,b1)).*exp(-b1);
        bd1=((n./b1-1).*besseli(n,b1)+besseli(n+1,b1)).*exp(-b1);
        bd2=((2+n.*(n-2.*b1-1)./b1.^2).*besseli(n,b1)-(2+1./b1).*besseli(n+1,b1)).*exp(-b1);
        %xx
        c1=n.^2./b1.*besseli(n,b1).*exp(-b1);
        c2=b1.*bd2+bd;
        SigLxx=SigLxx+Z(zn1).*(c1.*cos(de1).^2+c2.*sin(de1).^2);
        %xy
        c3=n.^2./b1.*besseli(n,b1).*exp(-b1)-b1.*bd2-bd;
        SigLxy=SigLxy-1i.*Z(zn1).*(n.*bd1+1i.*sin(de1).*cos(de1).*c3);
        %xz
        c4=1i.*n.*vepar.*besseli(n,b1).*exp(-b1);
        c5=veper.*bd1;
        SigLxz=SigLxz-2i.*(1+zn1.*Z(zn1)).*(c4.*cos(de1)+c5.*sin(de1));
        %yx
        SigLyx=SigLyx+1i.*Z(zn1).*(n.*bd1-1i.*sin(de1).*cos(de1).*c3);
        %yy
        SigLyy=SigLyy+Z(zn1).*(c1.*sin(de1).^2+c2.*cos(de1).^2);
        %yz
        SigLyz=SigLyz-2i.*(1+zn1.*Z(zn1)).*(c4.*sin(de1)-c5.*cos(de1));
        %zx
        SigLzx=SigLzx-2i.*(1+zn1.*Z(zn1)).*(c4.*cos(de1)-c5.*sin(de1));
        %zy
        SigLzy=SigLzy-2i.*(1+zn1.*Z(zn1)).*(c4.*sin(de1)+c5.*cos(de1));
        %zz
        SigLzz=SigLzz+2.*zn1.*besseli(n,b1).*exp(-b1).*(1+zn1.*Z(zn1));
    end
    y = wpe.^2./w1./n1z./ve.*[SigLxx SigLxy SigLxz
        SigLyx SigLyy SigLyz
        SigLzx SigLzy SigLzz];
    num = num+1;
end
y = y;