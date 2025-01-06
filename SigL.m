%电子 低频 线性
function y = SigL(w)
global wpe wce nz ve de nper;
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
    zetan=(w-n.*wce)./nz./ve;
    b=(nper.*ve./wce).^2./2;
    vepar = wce./nper./ve;
    veper=nper.*ve./2./wce;
    bd=(besseli(n+1,b)+n.*besseli(n,b)./b).*exp(-b);
    bd1=(besseli(n+1,b)+(n./b-1).*besseli(n,b)).*exp(-b);
    bd2=((2+n.*(n-2.*b-1)./b.^2).*besseli(n,b)-(2+1./b).*besseli(n+1,b)).*exp(-b);
    %xx
    c1=(n.^2.*besseli(n,b)./b).*exp(-b);
    c2=b.*bd2+bd;
    SigLxx=SigLxx+Z(zetan).*(c1.*cos(de).^2+c2.*sin(de).^2);
    %xy
    c3=(n.^2.*besseli(n,b)./b).*exp(-b)-b.*bd2-bd;
    SigLxy=SigLxy-1i.*Z(zetan).*(n.*bd1+1i.*sin(de).*cos(de).*c3);
    %xz
    c4=1i.*n.*vepar.*besseli(n,b).*exp(-b);
    c5=veper.*bd1;
    SigLxz=SigLxz-2i.*(1+zetan.*Z(zetan)).*(c4.*cos(de)+c5.*sin(de));
    %yx
    SigLyx=SigLyx+1i.*Z(zetan).*(n.*bd1-1i.*sin(de).*cos(de).*c3);
    %yy
    SigLyy=SigLyy+Z(zetan).*(c1.*sin(de).^2+c2.*cos(de).^2);
    %yz
    SigLyz=SigLyz-2i.*(1+zetan.*Z(zetan)).*(c4.*sin(de)-c5.*cos(de));
    %zx
    SigLzx=SigLzx-2i.*(1+zetan.*Z(zetan)).*(c4.*cos(de)-c5.*sin(de));
    %zy
    SigLzy=SigLzy-2i.*(1+zetan.*Z(zetan)).*(c4.*sin(de)+c5.*cos(de));
    %zz
    SigLzz=SigLzz+2.*zetan.*besseli(n,b).*exp(-b).*(1+zetan.*Z(zetan));
end
y = wpe.^2./w./nz./ve.*[SigLxx SigLxy SigLxz
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
        zetan=(w-n.*wce)./nz./ve;
        b=(nper.*ve./wce).^2./2;
        vepar = wce./nper./ve;
        veper=nper.*ve./2./wce;
        bd=(besseli(n+1,b)+n.*besseli(n,b)./b).*exp(-b);
        bd1=(besseli(n+1,b)+(n./b-1).*besseli(n,b)).*exp(-b);
        bd2=((2+n.*(n-2.*b-1)./b.^2).*besseli(n,b)-(2+1./b).*besseli(n+1,b)).*exp(-b);
        %xx
        c1=(n.^2.*besseli(n,b)./b).*exp(-b);
        c2=b.*bd2+bd;
        SigLxx=SigLxx+Z(zetan).*(c1.*cos(de).^2+c2.*sin(de).^2);
        %xy
        c3=(n.^2.*besseli(n,b)./b).*exp(-b)-b.*bd2-bd;
        SigLxy=SigLxy-1i.*Z(zetan).*(n.*bd1+1i.*sin(de).*cos(de).*c3);
        %xz
        c4=1i.*n.*vepar.*besseli(n,b).*exp(-b);
        c5=veper.*bd1;
        SigLxz=SigLxz-2i.*(1+zetan.*Z(zetan)).*(c4.*cos(de)+c5.*sin(de));
        %yx
        SigLyx=SigLyx+1i.*Z(zetan).*(n.*bd1-1i.*sin(de).*cos(de).*c3);
        %yy
        SigLyy=SigLyy+Z(zetan).*(c1.*sin(de).^2+c2.*cos(de).^2);
        %yz
        SigLyz=SigLyz-2i.*(1+zetan.*Z(zetan)).*(c4.*sin(de)-c5.*cos(de));
        %zx
        SigLzx=SigLzx-2i.*(1+zetan.*Z(zetan)).*(c4.*cos(de)-c5.*sin(de));
        %zy
        SigLzy=SigLzy-2i.*(1+zetan.*Z(zetan)).*(c4.*sin(de)+c5.*cos(de));
        %zz
        SigLzz=SigLzz+2.*zetan.*besseli(n,b).*exp(-b).*(1+zetan.*Z(zetan));
    end
    y = wpe.^2./w./nz./ve.*[SigLxx SigLxy SigLxz
        SigLyx SigLyy SigLyz
        SigLzx SigLzy SigLzz];
    num = num+1;
end
y = y;