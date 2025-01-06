%下边频 非线性
%完成于2015/9/10
function y = SigNL1(w)
w1=w-1;
global wpe wce ve;
global de de1 w0 nz n0z n1z nper n0x n1per nx ny;
global E0x E0y E0z;
E0xc=conj(E0x);
E0yc=conj(E0y);
E0zc=conj(E0z);
epsilon = 1e-09;
SigNLxx=0;
SigNLxy=0;
SigNLxz=0;
SigNLyx=0;
SigNLyy=0;
SigNLyz=0;
SigNLzx=0;
SigNLzy=0;
SigNLzz=0;
for n=-0:0
    for l=-0:0
        %等离子体色散函数
        zl=(w-l.*wce)./nz./ve;
        zn1=(w1-n.*wce)./n1z./ve;
        zl0=(1+(n-l).*wce)./n0z./ve;
        Z1=(Z(zl)-Z(zn1))./(zl-zn1);
        Z2=n1z.*Z(zl)-nz.*Z(zn1);
        Z3=(1+zn1.*Z(zn1))./(zl-zn1)+(Z(zl)-Z(zn1))./2./(zl-zn1).^2;
        Z4=(zl.*Z(zl)-zn1.*Z(zn1))./(zl-zn1);
        Z5=(Z(zl)-Z(zn1))./(zl-zn1).^2+2.*(1+zl.*Z(zl))./(zl-zn1);
        Z6=(Z(zl0)-Z(zn1))./(zl0-zn1);
        Z7=n1z.*Z(zl0)+n0z.*Z(zn1);
        Z8=-2.*(1+zl0.*Z(zl0))./(zl0-zn1)-(Z(zl0)-Z(zn1))./(zl0-zn1).^2;
        Z9=(zl0.*Z(zl0)-zn1.*Z(zn1))./(zl0-zn1);
        Z10=n1z.*(1+zl0.*Z(zl0))+n0z.*(1+zn1.*Z(zn1));
        Z11=1+(zl0.^2.*Z(zl0)-zn1.^2.*Z(zn1))./(zl0-zn1);
        Z12=n1z.*(1+zl.*Z(zl))-nz.*(1+zn1.*Z(zn1));
        Z13=1+(zl.^2.*Z(zl)-zn1.^2.*Z(zn1))./(zl-zn1);
        Z14=zn1.*Z3;
        Z15=2.*zl.*(1+zl.*Z(zl))./(zl-zn1)+zn1.*(Z(zl)-Z(zn1))./(zl-zn1).^2;
        Z16=2.*zl0.*(1+zl0.*Z(zl0))./(zl0-zn1)+zn1.*(Z(zl0)-Z(zn1))./(zl0-zn1).^2;
        Z17=(zl0.^2.*(1+zl0.*Z(zl0))-zn1.^2.*(1+zn1.*Z(zn1)))./(zl0-zn1);
        Z18=(zl.^2.*(1+zl.*Z(zl))-zn1.^2.*(1+zn1.*Z(zn1)))./(zl-zn1);
        %贝塞尔函数积分
        p=nper.*ve./wce;
        p0=n0x.*ve./wce;
        p1=n1per.*ve./wce;
        f1=@(x) x.*exp(-x.^2).*besselj(l,p.*x).*besselj(-(n-l),p0.*x).*besselj(n,p1.*x);
        i1=integral(f1,0,inf);
        f2=@(x) x.^2.*exp(-x.^2).*besselj(l,p.*x) ...
            .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)).*besselj(n,p1.*x);
        i2=integral(f2,0,inf);
        f3=@(x) x.^2.*exp(-x.^2).*besselj(l,p.*x).*besselj(-(n-l),p0.*x) ...
            .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
        i3=integral(f3,0,inf);
        f4=@(x) x.^3.*exp(-x.^2).*besselj(l,p.*x) ...
            .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)) ...
            .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
        i4=integral(f4,0,inf);
        f5=@(x) x.^2.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
            .*besselj(-(n-l),p0.*x).*besselj(n,p1.*x);
        i5=integral(f5,0,inf);
        f6=@(x) x.^3.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
            .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)).*besselj(n,p1.*x);
        i6=integral(f6,0,inf);
        f7=@(x) x.^3.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
            .*besselj(-(n-l),p0.*x).*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
        i7=integral(f7,0,inf);
        f8=@(x) x.^4.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
            .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)) ...
            .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
        i8=integral(f8,0,inf);
        f9=@(x) exp(-x.^2)./x.*besselj(l,p.*x).*besselj(-(n-l),p0.*x).*besselj(n,p1.*x);
        i9=integral(f9,1e-39,inf);
        f10=@(x) exp(-x.^2).*besselj(l,p.*x) ...
            .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)).*besselj(n,p1.*x);
        i10=integral(f10,0,inf);
        f11=@(x) exp(-x.^2).*besselj(l,p.*x).*besselj(-(n-l),p0.*x) ...
            .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
        i11=integral(f11,0,inf);
        f12=@(x) x.*exp(-x.^2).*besselj(l,p.*x) ...
            .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)) ...
            .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
        i12=integral(f12,0,inf);
        f13=@(x) exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
            .*besselj(-(n-l),p0.*x).*besselj(n,p1.*x);
        i13=integral(f13,0,inf);
        f14=@(x) x.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
            .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)).*besselj(n,p1.*x);
        i14=integral(f14,0,inf);
        f15=@(x) x.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
            .*besselj(-(n-l),p0.*x).*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
        i15=integral(f15,0,inf);
        f16=@(x) x.^2.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
            .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)) ...
            .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
        i16=integral(f16,0,inf);
        % 通用的项
        b1=-(n0x.*ny-2i.*(n-l).*wce.^2./ve.^2);
        b2=-(n0x.*ny+2i.*wce./ve.^2);
        b3=n0x.*ny+2i.*l.*wce.^2./ve.^2;
        b4=n0x.*ny+2i.*wce.*w./ve.^2;
% Part 1
        d1=-E0xc.*n.*(n-l).*wce.^2./n0x./n1per./ve.^2;
        d2=1i.*E0yc.*n.*wce./n1per./ve;
        d3=1i.*E0xc.*(n-l).*wce./n0x./ve;
        d4=E0yc;
        d7=-E0xc.*(n-l).*wce./n0x./ve;
        d8=1i.*E0yc;
        c1=b1.*l.*wce;
        c2=-(n-l).*nper.^2.*wce;
        c3=-2i.*l.*wce.^2./ve;
        c4=nper.^2.*ve;
        c5=2i.*l.*n0z.*wce.^2./ve;
        axx1=((c1.*cos(de)+c2.*sin(de)).*Z1+(c3.*cos(de)+c4.*sin(de)).*Z2+c5.*cos(de).*Z3)./nper ...
            .*(d1.*cos(de1).*i1+d2.*cos(de1).*i2+d3.*sin(de1).*i3+d4.*sin(de1).*i4);
        axy1=((c1.*sin(de)-c2.*cos(de)).*Z1+(c3.*sin(de)-c4.*cos(de)).*Z2+c5.*sin(de).*Z3)./nper ...
            .*(d1.*cos(de1).*i1+d2.*cos(de1).*i2+d3.*sin(de1).*i3+d4.*sin(de1).*i4);
        axz1=(b1.*Z4-2i.*wce./ve.*Z12+2i.*n0z.*wce./ve.*(-Z13+zl./2.*Z5)) ...
            .*(d1.*cos(de1).*i1+d2.*cos(de1).*i2+d3.*sin(de1).*i3+d4.*sin(de1).*i4);
        ayx1=((c1.*cos(de)+c2.*sin(de)).*Z1+(c3.*cos(de)+c4.*sin(de)).*Z2+c5.*cos(de).*Z3)./nper ...
            .*(d1.*sin(de1).*i1+d2.*sin(de1).*i2-d3.*cos(de1).*i3-d4.*cos(de1).*i4);
        ayy1=((c1.*sin(de)-c2.*cos(de)).*Z1+(c3.*sin(de)-c4.*cos(de)).*Z2+c5.*sin(de).*Z3)./nper ...
            .*(d1.*sin(de1).*i1+d2.*sin(de1).*i2-d3.*cos(de1).*i3-d4.*cos(de1).*i4);
        ayz1=(b1.*Z4-2i.*wce./ve.*Z12+2i.*n0z.*wce./ve.*(-Z13+zl./2.*Z5)) ...
            .*(d1.*sin(de1).*i1+d2.*sin(de1).*i2-d3.*cos(de1).*i3-d4.*cos(de1).*i4);
        azx1=((c1.*cos(de)+c2.*sin(de)).*Z4+(c3.*cos(de)+c4.*sin(de)).*Z12+c5.*cos(de).*Z14)./nper ...
            .*(d7.*i1+d8.*i2);
        azy1=((c1.*sin(de)-c2.*cos(de)).*Z4+(c3.*sin(de)-c4.*cos(de)).*Z12+c5.*sin(de).*Z14)./nper ...
            .*(d7.*i1+d8.*i2);
        azz1=(b2.*Z13+1i.*n0z.*wce./ve.*zl.*Z15).*(d7.*i1+d8.*i2);
        % Next
        c1=1i.*b1.*ve;
        c2=2.*wce;
        c3=-2.*n0z.*wce;
        axx2=sin(de).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*cos(de1).*i5+d2.*cos(de1).*i6+d3.*sin(de1).*i7+d4.*sin(de1).*i8);
        axy2=-cos(de).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*cos(de1).*i5+d2.*cos(de1).*i6+d3.*sin(de1).*i7+d4.*sin(de1).*i8);
        ayx2=sin(de).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*sin(de1).*i5+d2.*sin(de1).*i6-d3.*cos(de1).*i7-d4.*cos(de1).*i8);
        ayy2=-cos(de).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*sin(de1).*i5+d2.*sin(de1).*i6-d3.*cos(de1).*i7-d4.*cos(de1).*i8);
        azx2=sin(de).*(c1.*Z4+c2.*Z12+c3.*Z14).*(d7.*i5+d8.*i6);
        azy2=-cos(de).*(c1.*Z4+c2.*Z12+c3.*Z14).*(d7.*i5+d8.*i6);
        % Next
        c1=-l.^2.*wce.^2./nper./ve;
        c2=-(n-l).*wce./ve;
        axx3=c1.*sin(de).*(Z2+c2.*Z1).*(d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12);
        axy3=-c1.*cos(de).*(Z2+c2.*Z1).*(d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12);
        ayx3=c1.*sin(de).*(Z2+c2.*Z1).*(d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12);
        ayy3=-c1.*cos(de).*(Z2+c2.*Z1).*(d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12);
        azx3=c1.*sin(de).*(Z12+c2.*Z4).*(d7.*i9+d8.*i10);
        azy3=-c1.*cos(de).*(Z12+c2.*Z4).*(d7.*i9+d8.*i10);
        % Next
        c1=1i.*l.*wce;
        axx4=c1.*cos(de).*(Z2+c2.*Z1).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
        axy4=c1.*sin(de).*(Z2+c2.*Z1).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
        axz4=1i.*nper.*(Z12-(n-l).*wce./ve.*Z4).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
        ayx4=c1.*cos(de).*(Z2+c2.*Z1).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
        ayy4=c1.*sin(de).*(Z2+c2.*Z1).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
        ayz4=1i.*nper.*(Z12-(n-l).*wce./ve.*Z4).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
        azx4=c1.*cos(de).*(Z12+c2.*Z4).*(d7.*i13+d8.*i14);
        azy4=c1.*sin(de).*(Z12+c2.*Z4).*(d7.*i13+d8.*i14);
        azz4=1i.*nper.*(Z13./ve-n0z.*Z18).*(d7.*i13+d8.*i14);
        % Sum
        sxx1=30./w0./nz./n1z.*(axx1+axx2+axx3+axx4);
        sxy1=30./w0./nz./n1z.*(axy1+axy2+axy3+axy4);
        sxz1=30./w0./nz./n1z.*ve.*(axz1+axz4);
        syx1=30./w0./nz./n1z.*(ayx1+ayx2+ayx3+ayx4);
        syy1=30./w0./nz./n1z.*(ayy1+ayy2+ayy3+ayy4);
        syz1=30./w0./nz./n1z.*ve.*(ayz1+ayz4);
        szx1=30./w0./nz./n1z.*(azx1+azx2+azx3+azx4);
        szy1=30./w0./nz./n1z.*(azy1+azy2+azy3+azy4);
        szz1=30./w0./nz./n1z.*ve.*(azz1+azz4);
% Part 2
        d5=n.*wce./n1per./ve;
        c1=b2.*l.*wce;
        c2=-(n-l).*nper.^2.*wce;
        c3=1i.*l.*(1+(n-l).*wce).*wce.^2./ve.^2;
        axx1=((c1.*cos(de)+c2.*sin(de)).*Z4+c3.*cos(de).*Z5)./nper.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        axy1=((c1.*sin(de)-c2.*cos(de)).*Z4+c3.*sin(de).*Z5)./nper.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        axz1=(b2.*ve.*Z13+1i.*wce.*(1+(n-l).*wce)./ve.*zl.*Z5).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        ayx1=((c1.*cos(de)+c2.*sin(de)).*Z4+c3.*cos(de).*Z5)./nper.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        ayy1=((c1.*sin(de)-c2.*cos(de)).*Z4+c3.*sin(de).*Z5)./nper.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        ayz1=(b2.*ve.*Z13+1i.*wce.*(1+(n-l).*wce)./ve.*zl.*Z5).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        azx1=((c1.*cos(de)+c2.*sin(de)).*Z13+c3.*cos(de).*Z15)./nper.*i1;
        azy1=((c1.*sin(de)-c2.*cos(de)).*Z13+c3.*sin(de).*Z15)./nper.*i1;
        azz1=(b2.*ve.*Z18+1i.*wce.*(1+(n-l).*wce)./ve.*zl.*Z15).*i1;
        % Next
        c1=1i.*b2;
        c2=-(1+(n-l).*wce).*wce./ve.^2;
        axx2=ve.*sin(de).*(c1.*Z4+c2.*Z5).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
        axy2=-ve.*cos(de).*(c1.*Z4+c2.*Z5).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
        ayx2=ve.*sin(de).*(c1.*Z4+c2.*Z5).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
        ayy2=-ve.*cos(de).*(c1.*Z4+c2.*Z5).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
        azx2=ve.*sin(de).*(c1.*Z13+c2.*Z15).*i5;
        azy2=-ve.*cos(de).*(c1.*Z13+c2.*Z15).*i5;
        % Next
        c1=l.^2.*(n-l).*wce.^3./nper./ve.^2;
        axx3=c1.*sin(de).*Z4.*(d5.*cos(de1).*i9-1i.*sin(de1).*i11);
        axy3=-c1.*cos(de).*Z4.*(d5.*cos(de1).*i9-1i.*sin(de1).*i11);
        ayx3=c1.*sin(de).*Z4.*(d5.*sin(de1).*i9+1i.*cos(de1).*i11);
        ayy3=-c1.*cos(de).*Z4.*(d5.*sin(de1).*i9+1i.*cos(de1).*i11);
        azx3=c1.*sin(de).*Z13.*i9;
        azy3=-c1.*cos(de).*Z13.*i9;
        % Next
        c1=-1i.*l.*(n-l).*wce.^2./ve;
        axx4=c1.*cos(de).*Z4.*(d5.*cos(de1).*i13-1i.*sin(de1).*i15);
        axy4=c1.*sin(de).*Z4.*(d5.*cos(de1).*i13-1i.*sin(de1).*i15);
        axz4=-1i.*(n-l).*nper.*wce.*Z13.*(d5.*cos(de1).*i13-1i.*sin(de1).*i15);
        ayx4=c1.*cos(de).*Z4.*(d5.*sin(de1).*i13+1i.*cos(de1).*i15);
        ayy4=c1.*sin(de).*Z4.*(d5.*sin(de1).*i13+1i.*cos(de1).*i15);
        ayz4=-1i.*(n-l).*nper.*wce.*Z13.*(d5.*sin(de1).*i13+1i.*cos(de1).*i15);
        azx4=c1.*cos(de).*Z13.*i13;
        azy4=c1.*sin(de).*Z13.*i13;
        azz4=-1i.*(n-l).*nper.*wce.*Z18.*i13;
        % Sum
        sxx2=30./w0./nz./n1z.*E0zc.*(axx1+axx2+axx3+axx4);
        sxy2=30./w0./nz./n1z.*E0zc.*(axy1+axy2+axy3+axy4);
        sxz2=30./w0./nz./n1z.*E0zc.*(axz1+axz4);
        syx2=30./w0./nz./n1z.*E0zc.*(ayx1+ayx2+ayx3+ayx4);
        syy2=30./w0./nz./n1z.*E0zc.*(ayy1+ayy2+ayy3+ayy4);
        syz2=30./w0./nz./n1z.*E0zc.*(ayz1+ayz4);
        szx2=30./w0./nz./n1z.*E0zc.*(azx1+azx2+azx3+azx4);
        szy2=30./w0./nz./n1z.*E0zc.*(azy1+azy2+azy3+azy4);
        szz2=30./w0./nz./n1z.*E0zc.*(azz1+azz4);
% Part 3
        d6=l.*wce./nper./ve;
        axx1=d6.*cos(de).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        axy1=d6.*sin(de).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        ayx1=d6.*cos(de).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        ayy1=d6.*sin(de).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        axx2=1i.*sin(de).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
        axy2=-1i.*cos(de).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
        ayx2=1i.*sin(de).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
        ayy2=-1i.*cos(de).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
        % Sum
        sxx3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z2.*(axx1+axx2);
        sxy3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z2.*(axy1+axy2);
        sxz3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z12.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        syx3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z2.*(ayx1+ayx2);
        syy3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z2.*(ayy1+ayy2);
        syz3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z12.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        szx3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z12.*(d6.*cos(de).*i1+1i.*sin(de).*i5);
        szy3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z12.*(d6.*sin(de).*i1-1i.*cos(de).*i5);
        szz3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*n0z.*ve.*(zl0.*Z13-Z18).*i1;
        % Next
        d1=-1i.*E0yc.*n.*(n-l).*wce.^2./n0x./n1per./ve.^2;
        d2=E0xc.*n.*wce./n1per./ve;
        d3=-E0yc.*(n-l).*wce./n0x./ve;
        d4=-1i.*E0xc;
        d7=-1i.*E0yc.*(n-l).*wce./n0x./ve;
        axx1=d6.*cos(de).*(d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12);
        axy1=d6.*sin(de).*(d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12);
        axz1=d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12;
        ayx1=d6.*cos(de).*(d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12);
        ayy1=d6.*sin(de).*(d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12);
        ayz1=d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12;
        azx1=d6.*cos(de).*(E0xc.*i10+d7.*i9);
        azy1=d6.*sin(de).*(E0xc.*i10+d7.*i9);
        azz1=E0xc.*i10+d7.*i9;
        axx2=1i.*sin(de).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
        axy2=-1i.*cos(de).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
        ayx2=1i.*sin(de).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
        ayy2=-1i.*cos(de).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
        azx2=1i.*sin(de).*(E0xc.*i14+d7.*i13);
        azy2=-1i.*cos(de).*(E0xc.*i14+d7.*i13);
        % Sum
        sxx4=30./w0./nz./n1z.*1i.*l.*wce.*(Z2-(n-l).*wce./ve.*Z1).*(axx1+axx2);
        sxy4=30./w0./nz./n1z.*1i.*l.*wce.*(Z2-(n-l).*wce./ve.*Z1).*(axy1+axy2);
        sxz4=30./w0./nz./n1z.*1i.*l.*wce.*(Z12-(n-l).*wce./ve.*Z4).*axz1;
        syx4=30./w0./nz./n1z.*1i.*l.*wce.*(Z2-(n-l).*wce./ve.*Z1).*(ayx1+ayx2);
        syy4=30./w0./nz./n1z.*1i.*l.*wce.*(Z2-(n-l).*wce./ve.*Z1).*(ayy1+ayy2);
        syz4=30./w0./nz./n1z.*1i.*l.*wce.*(Z12-(n-l).*wce./ve.*Z4).*ayz1;
        szx4=30./w0./nz./n1z.*1i.*l.*wce.*(Z12-(n-l).*wce./ve.*Z4).*(azx1+azx2);
        szy4=30./w0./nz./n1z.*1i.*l.*wce.*(Z12-(n-l).*wce./ve.*Z4).*(azy1+azy2);
        szz4=30./w0./nz./n1z.*1i.*l.*wce.*(Z13./ve-n0z.*Z18).*azz1;
        % Next
        axx1=d6.*cos(de).*(d5.*cos(de1).*i10-1i.*sin(de1).*i12);
        axy1=d6.*sin(de).*(d5.*cos(de1).*i10-1i.*sin(de1).*i12);
        ayx1=d6.*cos(de).*(d5.*sin(de1).*i10+1i.*cos(de1).*i12);
        ayy1=d6.*sin(de).*(d5.*sin(de1).*i10+1i.*cos(de1).*i12);
        axx2=1i.*sin(de).*(d5.*cos(de1).*i14-1i.*sin(de1).*i16);
        axy2=-1i.*cos(de).*(d5.*cos(de1).*i14-1i.*sin(de1).*i16);
        ayx2=1i.*sin(de).*(d5.*sin(de1).*i14+1i.*cos(de1).*i16);
        ayy2=-1i.*cos(de).*(d5.*sin(de1).*i14+1i.*cos(de1).*i16);
        % Sum
        sxx5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z4.*(axx1+axx2);
        sxy5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z4.*(axy1+axy2);
        sxz5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z13.*(d5.*cos(de1).*i10-1i.*sin(de1).*i12);
        syx5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z4.*(ayx1+ayx2);
        syy5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z4.*(ayy1+ayy2);
        syz5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z13.*(d5.*sin(de1).*i10+1i.*cos(de1).*i12);
        szx5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z13.*(d6.*cos(de).*i10+1i.*sin(de).*i14);
        szy5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z13.*(d6.*sin(de).*i10-1i.*cos(de).*i14);
        szz5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z18.*i10;
        % Next
        axx1=d6.*cos(de).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        axy1=d6.*sin(de).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        ayx1=d6.*cos(de).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        ayy1=d6.*sin(de).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        axx2=1i.*sin(de).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
        axy2=-1i.*cos(de).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
        ayx2=1i.*sin(de).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
        ayy2=-1i.*cos(de).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
        % Sum
        sxx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z1.*(axx1+axx2);
        sxy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z1.*(axy1+axy2);
        sxz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z4.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        syx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z1.*(ayx1+ayx2);
        syy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z1.*(ayy1+ayy2);
        syz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z4.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        szx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z4.*(d6.*cos(de).*i1+1i.*sin(de).*i5);
        szy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z4.*(d6.*sin(de).*i1-1i.*cos(de).*i5);
        szz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z13.*i1;
% Part 4
        d1=n.*l.*wce.^2./nper./n1per./ve.^2;
        d2=1i.*n.*wce./n1per./ve;
        d3=-1i.*l.*wce./nper./ve;
        c1=b3.*E0zc-2i.*(n-l).*nz.*wce.^2./n0x./ve.^2.*E0xc;
        c2=(-(n-l).*b3./n0x.*E0xc-l.*n0x.*E0yc).*wce./ve;
        c3=-(-2i.*(n-l).*wce.^2./n0x./ve.^2.*E0xc-n0x.*E0yc);
        c4=-2i.*wce./ve.*E0zc;
        c5=2i.*nz.*wce./ve.*E0zc;
        c6=1i.*wce.*((1+(n-l).*wce).*E0zc-(n-l).*n0z.*wce./n0x.*E0xc).*nz./n0z./ve.^2;
        c7=c1.*Z9+c2.*Z6+c3.*Z7+c4.*Z10+c5.*Z11+c6.*Z8;
        c8=b4.*ve;
        c9=-l.*n0x.*wce.*E0yc;
        c10=1i.*wce.*(w-l.*wce).*((1+(n-l).*wce).*E0zc-(n-l).*n0z.*wce./n0x.*E0xc)./n0z./ve.^2;
        c11=c8.*(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9)+c9.*Z9+c10.*Z8;
        c12=(n0x.*ny+2i.*w.*wce./ve.^2).*E0zc-2i.*(n-l).*nz.*wce.^2./n0x./ve.^2.*E0xc;
        c13=(-b4.*(n-l)./n0x.*E0xc-l.*n0x.*E0yc).*wce./ve;
        c14=-1i.*wce.*((1+(n-l).*wce).*E0zc-(n-l).*n0z.*wce./n0x.*E0xc).*nz./n0z./ve.^2;
        axx1=c7.*(d1.*cos(de).*cos(de1).*i1+d2.*sin(de).*cos(de1).*i5+d3.*cos(de).*sin(de1).*i3+sin(de).*sin(de1).*i7);
        axy1=c7.*(d1.*sin(de).*cos(de1).*i1-d2.*cos(de).*cos(de1).*i5+d3.*sin(de).*sin(de1).*i3-cos(de).*sin(de1).*i7);
        axz1=c11.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        ayx1=c7.*(d1.*cos(de).*sin(de1).*i1+d2.*sin(de).*sin(de1).*i5-d3.*cos(de).*cos(de1).*i3-sin(de).*cos(de1).*i7);
        ayy1=c7.*(d1.*sin(de).*sin(de1).*i1-d2.*cos(de).*sin(de1).*i5-d3.*sin(de).*cos(de1).*i3+cos(de).*cos(de1).*i7);
        ayz1=c11.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        azx1=(c12.*Z11+c2.*Z9+c3.*Z10-c6.*Z16).*(d6.*cos(de).*i1+1i.*sin(de).*i5);
        azy1=(c12.*Z11+c2.*Z9+c3.*Z10-c6.*Z16).*(d6.*sin(de).*i1-1i.*cos(de).*i5);
        azz1=(b4.*E0zc.*Z17+c13.*Z11+c14.*zl.*Z16).*i1;
        % Next
        c1=1i.*b3;
        c2=2.*wce./ve;
        c3=-nz.*wce./ve;
        c4=E0yc.*(c1.*Z6+c2.*Z7+c3.*(2.*Z9+Z8));
        axx2=c4.*(d1.*cos(de).*cos(de1).*i2+d2.*sin(de).*cos(de1).*i6+d3.*cos(de).*sin(de1).*i4+sin(de).*sin(de1).*i8);
        axy2=c4.*(d1.*sin(de).*cos(de1).*i2-d2.*cos(de).*cos(de1).*i6+d3.*sin(de).*sin(de1).*i4-cos(de).*sin(de1).*i8);
        axz2=ve.*E0yc.*(1i.*b4.*Z9-wce.*(w-l.*wce)./ve.^2.*Z8).*(d5.*cos(de1).*i2-1i.*sin(de1).*i4);
        ayx2=c4.*(d1.*cos(de).*sin(de1).*i2+d2.*sin(de).*sin(de1).*i6-d3.*cos(de).*cos(de1).*i4-sin(de).*cos(de1).*i8);
        ayy2=c4.*(d1.*sin(de).*sin(de1).*i2-d2.*cos(de).*sin(de1).*i6-d3.*sin(de).*cos(de1).*i4+cos(de).*cos(de1).*i8);
        ayz2=ve.*E0yc.*(1i.*b4.*Z9-wce.*(w-l.*wce)./ve.^2.*Z8).*(d5.*sin(de1).*i2+1i.*cos(de1).*i4);
        azx2=E0yc.*(c1.*Z9+c2.*Z10+c3.*(2.*Z11-Z16)).*(d6.*cos(de).*i2+1i.*sin(de).*i6);
        azy2=E0yc.*(c1.*Z9+c2.*Z10+c3.*(2.*Z11-Z16)).*(d6.*sin(de).*i2-1i.*cos(de).*i6);
        azz2=E0yc.*(1i.*b4.*Z11+nz.*wce./ve.*zl.*Z16).*i2;
        % Next
        c1=-(n-l).^2.*wce.^2./n0x./ve.^2.*E0yc.*(Z7-l.*wce./ve.*Z6);
        axx3=c1.*(d1.*cos(de).*cos(de1).*i9+d2.*sin(de).*cos(de1).*i13+d3.*cos(de).*sin(de1).*i11+sin(de).*sin(de1).*i15);
        axy3=c1.*(d1.*sin(de).*cos(de1).*i9-d2.*cos(de).*cos(de1).*i13+d3.*sin(de).*sin(de1).*i11-cos(de).*sin(de1).*i15);
        axz3=E0yc.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z9.*(d5.*cos(de1).*i9-1i.*sin(de1).*i11);
        ayx3=c1.*(d1.*cos(de).*sin(de1).*i9+d2.*sin(de).*sin(de1).*i13-d3.*cos(de).*cos(de1).*i11-sin(de).*cos(de1).*i15);
        ayy3=c1.*(d1.*sin(de).*sin(de1).*i9-d2.*cos(de).*sin(de1).*i13-d3.*sin(de).*cos(de1).*i11+cos(de).*cos(de1).*i15);
        ayz3=E0yc.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z9.*(d5.*sin(de1).*i9+1i.*cos(de1).*i11);
        azx3=-(n-l).^2.*wce.^2./n0x./ve.^2.*E0yc.*(Z10-l.*wce./ve.*Z9).*(d6.*cos(de).*i9+1i.*sin(de).*i13);
        azy3=-(n-l).^2.*wce.^2./n0x./ve.^2.*E0yc.*(Z10-l.*wce./ve.*Z9).*(d6.*sin(de).*i9-1i.*cos(de).*i13);
        azz3=l.*(n-l).^2.*wce.^3./n0x./ve.^3.*E0yc.*Z11.*i9;
        % Next
        c1=1i.*n0x.*E0zc;
        c2=-1i.*l.*n0x.*wce./ve.*E0zc;
        c3=-1i.*(n-l).*wce./ve.*E0xc;
        c4=1i.*l.*(n-l).*wce.^2./ve.^2.*E0xc;
        c5=c1.*Z10+c2.*Z9+c3.*Z7+c4.*Z6;
        c6=-1i.*n0x.*w./ve.*E0zc;
        c7=1i.*nz.*n0x.*E0zc;
        axx4=c5.*(d1.*cos(de).*cos(de1).*i10+d2.*sin(de).*cos(de1).*i14+d3.*cos(de).*sin(de1).*i12+sin(de).*sin(de1).*i16);
        axy4=c5.*(d1.*sin(de).*cos(de1).*i10-d2.*cos(de).*cos(de1).*i14+d3.*sin(de).*sin(de1).*i12-cos(de).*sin(de1).*i16);
        axz4=-1i.*l.*n0x.*wce.*(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*(d5.*cos(de1).*i10-1i.*sin(de1).*i12);
        ayx4=c5.*(d1.*cos(de).*sin(de1).*i10+d2.*sin(de).*sin(de1).*i14-d3.*cos(de).*cos(de1).*i12-sin(de).*cos(de1).*i16);
        ayy4=c5.*(d1.*sin(de).*sin(de1).*i10-d2.*cos(de).*sin(de1).*i14-d3.*sin(de).*cos(de1).*i12+cos(de).*cos(de1).*i16);
        ayz4=-1i.*l.*n0x.*wce.*(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*(d5.*sin(de1).*i10+1i.*cos(de1).*i12);
        azx4=(c6.*Z11+c7.*Z17+c3.*Z10+c4.*Z9).*(d6.*cos(de).*i10+1i.*sin(de).*i14);
        azy4=(c6.*Z11+c7.*Z17+c3.*Z10+c4.*Z9).*(d6.*sin(de).*i10-1i.*cos(de).*i14);
        azz4=-1i.*l.*n0x.*wce./ve.*(E0zc.*Z17-E0xc.*(n-l).*wce./n0x./ve.*Z11).*i10;
        % Sum
        sxx7=30./w0./n0z./n1z.*ve./w.*(axx1+axx2+axx3+axx4);
        sxy7=30./w0./n0z./n1z.*ve./w.*(axy1+axy2+axy3+axy4);
        sxz7=30./w0./n0z./n1z./w.*(axz1+axz2+axz3+axz4);
        syx7=30./w0./n0z./n1z.*ve./w.*(ayx1+ayx2+ayx3+ayx4);
        syy7=30./w0./n0z./n1z.*ve./w.*(ayy1+ayy2+ayy3+ayy4);
        syz7=30./w0./n0z./n1z./w.*(ayz1+ayz2+ayz3+ayz4);
        szx7=30./w0./n0z./n1z.*ve./w.*(azx1+azx2+azx3+azx4);
        szy7=30./w0./n0z./n1z.*ve./w.*(azy1+azy2+azy3+azy4);
        szz7=30./w0./n0z./n1z.*ve./w.*(azz1+azz2+azz3+azz4);
% Part 5
        axy1=(E0zc.*Z10-E0xc.*(n-l).*wce./n0x./ve.*Z7).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        ayy1=(E0zc.*Z10-E0xc.*(n-l).*wce./n0x./ve.*Z7).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        azy1=(E0zc.*((w-l.*wce)./ve.*Z11-nz.*Z17)-E0xc.*(n-l).*wce./n0x./ve.*(-Z10)).*i1;
        axy2=1i.*E0yc.*Z7.*(d5.*cos(de1).*i2-1i.*sin(de1).*i4);
        ayy2=1i.*E0yc.*Z7.*(d5.*sin(de1).*i2+1i.*cos(de1).*i4);
        azy2=-1i.*E0yc.*Z10.*i2;
        % Sum
        sxy8=-30./w0./n0z./n1z.*n0x.*ve./w.*(axy1+axy2);
        syy8=-30./w0./n0z./n1z.*n0x.*ve./w.*(ayy1+ayy2);
        szy8=30./w0./n0z./n1z.*n0x.*ve./w.*(azy1+azy2);
        % Next
        c1=E0zc.*(-Z10+l.*wce./ve.*Z9)-E0xc.*(n-l).*wce./n0x./ve.*(-Z7+l.*wce./ve.*Z6);
        axx1=c1.*(-d1.*sin(de).*cos(de1).*i9+d2.*cos(de).*cos(de1).*i13-d3.*sin(de).*sin(de1).*i11+cos(de).*sin(de1).*i15);
        axy1=c1.*(d1.*cos(de).*cos(de1).*i9+d2.*sin(de).*cos(de1).*i13+d3.*cos(de).*sin(de1).*i11+sin(de).*sin(de1).*i15);
        axz1=(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*(d5.*cos(de1).*i13-1i.*sin(de1).*i15);
        ayx1=c1.*(-d1.*sin(de).*sin(de1).*i9+d2.*cos(de).*sin(de1).*i13+d3.*sin(de).*cos(de1).*i11-cos(de).*cos(de1).*i15);
        ayy1=c1.*(d1.*cos(de).*sin(de1).*i9+d2.*sin(de).*sin(de1).*i13-d3.*cos(de).*cos(de1).*i11-sin(de).*cos(de1).*i15);
        ayz1=(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*(d5.*sin(de1).*i13+1i.*cos(de1).*i15);
        azx1=(E0zc.*(w./ve.*Z11-nz.*Z17)-E0xc.*(n-l).*wce./n0x./ve.*(-Z10+l.*wce./ve.*Z9)).*(d6.*sin(de).*i9-1i.*cos(de).*i13);
        azy1=(E0zc.*(w./ve.*Z11-nz.*Z17)-E0xc.*(n-l).*wce./n0x./ve.*(-Z10+l.*wce./ve.*Z9)).*(-d6.*cos(de).*i9-1i.*sin(de).*i13);
        azz1=(E0zc.*Z17-E0xc.*(n-l).*wce./n0x./ve.*Z11).*i13;
        c1=1i.*E0yc.*(-Z7+l.*wce./ve.*Z6);
        axx2=c1.*(-d1.*sin(de).*cos(de1).*i10+d2.*cos(de).*cos(de1).*i14-d3.*sin(de).*sin(de1).*i12+cos(de).*sin(de1).*i16);
        axy2=c1.*(d1.*cos(de).*cos(de1).*i10+d2.*sin(de).*cos(de1).*i14+d3.*cos(de).*sin(de1).*i12+sin(de).*sin(de1).*i16);
        axz2=1i.*E0yc.*Z9.*(d5.*cos(de1).*i14-1i.*sin(de1).*i16);
        ayx2=c1.*(-d1.*sin(de).*sin(de1).*i10+d2.*cos(de).*sin(de1).*i14+d3.*sin(de).*cos(de1).*i12-cos(de).*cos(de1).*i16);
        ayy2=c1.*(d1.*cos(de).*sin(de1).*i10+d2.*sin(de).*sin(de1).*i14-d3.*cos(de).*cos(de1).*i12-sin(de).*cos(de1).*i16);
        ayz2=1i.*E0yc.*Z9.*(d5.*sin(de1).*i14+1i.*cos(de1).*i16);
        azx2=1i.*E0yc.*(-Z10+l.*wce./ve.*Z9).*(d6.*sin(de).*i10-1i.*cos(de).*i14);
        azy2=1i.*E0yc.*(-Z10+l.*wce./ve.*Z9).*(-d6.*cos(de).*i10-1i.*sin(de).*i14);
        azz2=1i.*E0yc.*Z11.*i14;
        % Sum
        sxx9=30./w0./n0z./n1z.*(n-l).*wce./w.*(axx1+axx2);
        sxy9=30./w0./n0z./n1z.*(n-l).*wce./w.*(axy1+axy2);
        sxz9=30./w0./n0z./n1z.*1i.*(n-l).*nper.*wce./w.*(axz1+axz2);
        syx9=30./w0./n0z./n1z.*(n-l).*wce./w.*(ayx1+ayx2);
        syy9=30./w0./n0z./n1z.*(n-l).*wce./w.*(ayy1+ayy2);
        syz9=30./w0./n0z./n1z.*1i.*(n-l).*nper.*wce./w.*(ayz1+ayz2);
        szx9=-30./w0./n0z./n1z.*(n-l).*wce./w.*(azx1+azx2);
        szy9=-30./w0./n0z./n1z.*(n-l).*wce./w.*(azy1+azy2);
        szz9=30./w0./n0z./n1z.*1i.*(n-l).*nper.*wce./w.*(azz1+azz2);
        % Next
        ax1=(E0zc.*Z9-E0xc.*(n-l).*wce./n0x./ve.*Z6).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
        ay1=(E0zc.*Z9-E0xc.*(n-l).*wce./n0x./ve.*Z6).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
        az1=(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*i1;
        ax2=1i.*E0yc.*Z6.*(d5.*cos(de1).*i2-1i.*sin(de1).*i4);
        ay2=1i.*E0yc.*Z6.*(d5.*sin(de1).*i2+1i.*cos(de1).*i4);
        az2=1i.*E0yc.*Z9.*i2;
        % Sum
        sxx10=30./w0./n0z./n1z.*(n-l).*ny.*wce./w.*(ax1+ax2);
        sxy10=-30./w0./n0z./n1z.*(n-l).*nx.*wce./w.*(ax1+ax2);
        syx10=30./w0./n0z./n1z.*(n-l).*ny.*wce./w.*(ay1+ay2);
        syy10=-30./w0./n0z./n1z.*(n-l).*nx.*wce./w.*(ay1+ay2);
        szx10=30./w0./n0z./n1z.*(n-l).*ny.*wce./w.*(az1+az2);
        szy10=-30./w0./n0z./n1z.*(n-l).*nx.*wce./w.*(az1+az2);
        % % %
        sxx=sxx1+sxx2+sxx3+sxx4+sxx5+sxx6+sxx7+sxx9+sxx10;
        sxy=sxy1+sxy2+sxy3+sxy4+sxy5+sxy6+sxy7+sxy8+sxy9+sxy10;
        sxz=sxz1+sxz2+sxz3+sxz4+sxz5+sxz6+sxz7+sxz9;
        syx=syx1+syx2+syx3+syx4+syx5+syx6+syx7+syx9+syx10;
        syy=syy1+syy2+syy3+syy4+syy5+syy6+syy7+syy8+syy9+syy10;
        syz=syz1+syz2+syz3+syz4+syz5+syz6+syz7+syz9;
        szx=szx1+szx2+szx3+szx4+szx5+szx6+szx7+szx9+szx10;
        szy=szy1+szy2+szy3+szy4+szy5+szy6+szy7+szy8+szy9+szy10;
        szz=szz1+szz2+szz3+szz4+szz5+szz6+szz7+szz9;
        SigNLxx=SigNLxx+exp(1i.*(n.*de1-l.*de)).*sxx;
        SigNLxy=SigNLxy+exp(1i.*(n.*de1-l.*de)).*sxy;
        SigNLxz=SigNLxz+exp(1i.*(n.*de1-l.*de)).*sxz;
        SigNLyx=SigNLyx+exp(1i.*(n.*de1-l.*de)).*syx;
        SigNLyy=SigNLyy+exp(1i.*(n.*de1-l.*de)).*syy;
        SigNLyz=SigNLyz+exp(1i.*(n.*de1-l.*de)).*syz;
        SigNLzx=SigNLzx+exp(1i.*(n.*de1-l.*de)).*szx;
        SigNLzy=SigNLzy+exp(1i.*(n.*de1-l.*de)).*szy;
        SigNLzz=SigNLzz+exp(1i.*(n.*de1-l.*de)).*szz;
    end
end
%{
y=[SigNLxx SigNLxy SigNLxz
    SigNLyx SigNLyy SigNLyz
    SigNLzx SigNLzy SigNLzz];
%}
y = 0.096./81.9.*wpe.^2./w1./wce./ve.^2.*[SigNLxx SigNLxy SigNLxz
    SigNLyx SigNLyy SigNLyz
    SigNLzx SigNLzy SigNLzz];
yold = y-100;
num = 1;
while max(max(abs(y-yold))) >= epsilon
    yold = y;
    SigNLxx=0;
    SigNLxy=0;
    SigNLxz=0;
    SigNLyx=0;
    SigNLyy=0;
    SigNLyz=0;
    SigNLzx=0;
    SigNLzy=0;
    SigNLzz=0;
    for n=-num:num
        for l=-num:num
            %等离子体色散函数
            zl=(w-l.*wce)./nz./ve;
            zn1=(w1-n.*wce)./n1z./ve;
            zl0=(1+(n-l).*wce)./n0z./ve;
            Z1=(Z(zl)-Z(zn1))./(zl-zn1);
            Z2=n1z.*Z(zl)-nz.*Z(zn1);
            Z3=(1+zn1.*Z(zn1))./(zl-zn1)+(Z(zl)-Z(zn1))./2./(zl-zn1).^2;
            Z4=(zl.*Z(zl)-zn1.*Z(zn1))./(zl-zn1);
            Z5=(Z(zl)-Z(zn1))./(zl-zn1).^2+2.*(1+zl.*Z(zl))./(zl-zn1);
            Z6=(Z(zl0)-Z(zn1))./(zl0-zn1);
            Z7=n1z.*Z(zl0)+n0z.*Z(zn1);
            Z8=-2.*(1+zl0.*Z(zl0))./(zl0-zn1)-(Z(zl0)-Z(zn1))./(zl0-zn1).^2;
            Z9=(zl0.*Z(zl0)-zn1.*Z(zn1))./(zl0-zn1);
            Z10=n1z.*(1+zl0.*Z(zl0))+n0z.*(1+zn1.*Z(zn1));
            Z11=1+(zl0.^2.*Z(zl0)-zn1.^2.*Z(zn1))./(zl0-zn1);
            Z12=n1z.*(1+zl.*Z(zl))-nz.*(1+zn1.*Z(zn1));
            Z13=1+(zl.^2.*Z(zl)-zn1.^2.*Z(zn1))./(zl-zn1);
            Z14=zn1.*Z3;
            Z15=2.*zl.*(1+zl.*Z(zl))./(zl-zn1)+zn1.*(Z(zl)-Z(zn1))./(zl-zn1).^2;
            Z16=2.*zl0.*(1+zl0.*Z(zl0))./(zl0-zn1)+zn1.*(Z(zl0)-Z(zn1))./(zl0-zn1).^2;
            Z17=(zl0.^2.*(1+zl0.*Z(zl0))-zn1.^2.*(1+zn1.*Z(zn1)))./(zl0-zn1);
            Z18=(zl.^2.*(1+zl.*Z(zl))-zn1.^2.*(1+zn1.*Z(zn1)))./(zl-zn1);
            %贝塞尔函数积分
            p=nper.*ve./wce;
            p0=n0x.*ve./wce;
            p1=n1per.*ve./wce;
            f1=@(x) x.*exp(-x.^2).*besselj(l,p.*x).*besselj(-(n-l),p0.*x).*besselj(n,p1.*x);
            i1=integral(f1,0,inf);
            f2=@(x) x.^2.*exp(-x.^2).*besselj(l,p.*x) ...
                .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)).*besselj(n,p1.*x);
            i2=integral(f2,0,inf);
            f3=@(x) x.^2.*exp(-x.^2).*besselj(l,p.*x).*besselj(-(n-l),p0.*x) ...
                .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
            i3=integral(f3,0,inf);
            f4=@(x) x.^3.*exp(-x.^2).*besselj(l,p.*x) ...
                .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)) ...
                .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
            i4=integral(f4,0,inf);
            f5=@(x) x.^2.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
                .*besselj(-(n-l),p0.*x).*besselj(n,p1.*x);
            i5=integral(f5,0,inf);
            f6=@(x) x.^3.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
                .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)).*besselj(n,p1.*x);
            i6=integral(f6,0,inf);
            f7=@(x) x.^3.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
                .*besselj(-(n-l),p0.*x).*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
            i7=integral(f7,0,inf);
            f8=@(x) x.^4.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
                .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)) ...
                .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
            i8=integral(f8,0,inf);
            f9=@(x) exp(-x.^2)./x.*besselj(l,p.*x).*besselj(-(n-l),p0.*x).*besselj(n,p1.*x);
            i9=integral(f9,1e-39,inf);
            f10=@(x) exp(-x.^2).*besselj(l,p.*x) ...
                .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)).*besselj(n,p1.*x);
            i10=integral(f10,0,inf);
            f11=@(x) exp(-x.^2).*besselj(l,p.*x).*besselj(-(n-l),p0.*x) ...
                .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
            i11=integral(f11,0,inf);
            f12=@(x) x.*exp(-x.^2).*besselj(l,p.*x) ...
                .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)) ...
                .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
            i12=integral(f12,0,inf);
            f13=@(x) exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
                .*besselj(-(n-l),p0.*x).*besselj(n,p1.*x);
            i13=integral(f13,0,inf);
            f14=@(x) x.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
                .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)).*besselj(n,p1.*x);
            i14=integral(f14,0,inf);
            f15=@(x) x.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
                .*besselj(-(n-l),p0.*x).*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
            i15=integral(f15,0,inf);
            f16=@(x) x.^2.*exp(-x.^2).*(l.*besselj(l,p.*x)./(p.*x)-besselj(l+1,p.*x)) ...
                .*(-(n-l).*besselj(-(n-l),p0.*x)./(p0.*x)-besselj(-(n-l)+1,p0.*x)) ...
                .*(n.*besselj(n,p1.*x)./(p1.*x)-besselj(n+1,p1.*x));
            i16=integral(f16,0,inf);
            % 通用的项
            b1=-(n0x.*ny-2i.*(n-l).*wce.^2./ve.^2);
            b2=-(n0x.*ny+2i.*wce./ve.^2);
            b3=n0x.*ny+2i.*l.*wce.^2./ve.^2;
            b4=n0x.*ny+2i.*wce.*w./ve.^2;
    % Part 1
            d1=-E0xc.*n.*(n-l).*wce.^2./n0x./n1per./ve.^2;
            d2=1i.*E0yc.*n.*wce./n1per./ve;
            d3=1i.*E0xc.*(n-l).*wce./n0x./ve;
            d4=E0yc;
            d7=-E0xc.*(n-l).*wce./n0x./ve;
            d8=1i.*E0yc;
            c1=b1.*l.*wce;
            c2=-(n-l).*nper.^2.*wce;
            c3=-2i.*l.*wce.^2./ve;
            c4=nper.^2.*ve;
            c5=2i.*l.*n0z.*wce.^2./ve;
            axx1=((c1.*cos(de)+c2.*sin(de)).*Z1+(c3.*cos(de)+c4.*sin(de)).*Z2+c5.*cos(de).*Z3)./nper ...
                .*(d1.*cos(de1).*i1+d2.*cos(de1).*i2+d3.*sin(de1).*i3+d4.*sin(de1).*i4);
            axy1=((c1.*sin(de)-c2.*cos(de)).*Z1+(c3.*sin(de)-c4.*cos(de)).*Z2+c5.*sin(de).*Z3)./nper ...
                .*(d1.*cos(de1).*i1+d2.*cos(de1).*i2+d3.*sin(de1).*i3+d4.*sin(de1).*i4);
            axz1=(b1.*Z4-2i.*wce./ve.*Z12+2i.*n0z.*wce./ve.*(-Z13+zl./2.*Z5)) ...
                .*(d1.*cos(de1).*i1+d2.*cos(de1).*i2+d3.*sin(de1).*i3+d4.*sin(de1).*i4);
            ayx1=((c1.*cos(de)+c2.*sin(de)).*Z1+(c3.*cos(de)+c4.*sin(de)).*Z2+c5.*cos(de).*Z3)./nper ...
                .*(d1.*sin(de1).*i1+d2.*sin(de1).*i2-d3.*cos(de1).*i3-d4.*cos(de1).*i4);
            ayy1=((c1.*sin(de)-c2.*cos(de)).*Z1+(c3.*sin(de)-c4.*cos(de)).*Z2+c5.*sin(de).*Z3)./nper ...
                .*(d1.*sin(de1).*i1+d2.*sin(de1).*i2-d3.*cos(de1).*i3-d4.*cos(de1).*i4);
            ayz1=(b1.*Z4-2i.*wce./ve.*Z12+2i.*n0z.*wce./ve.*(-Z13+zl./2.*Z5)) ...
                .*(d1.*sin(de1).*i1+d2.*sin(de1).*i2-d3.*cos(de1).*i3-d4.*cos(de1).*i4);
            azx1=((c1.*cos(de)+c2.*sin(de)).*Z4+(c3.*cos(de)+c4.*sin(de)).*Z12+c5.*cos(de).*Z14)./nper ...
                .*(d7.*i1+d8.*i2);
            azy1=((c1.*sin(de)-c2.*cos(de)).*Z4+(c3.*sin(de)-c4.*cos(de)).*Z12+c5.*sin(de).*Z14)./nper ...
                .*(d7.*i1+d8.*i2);
            azz1=(b2.*Z13+1i.*n0z.*wce./ve.*zl.*Z15).*(d7.*i1+d8.*i2);
            % Next
            c1=1i.*b1.*ve;
            c2=2.*wce;
            c3=-2.*n0z.*wce;
            axx2=sin(de).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*cos(de1).*i5+d2.*cos(de1).*i6+d3.*sin(de1).*i7+d4.*sin(de1).*i8);
            axy2=-cos(de).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*cos(de1).*i5+d2.*cos(de1).*i6+d3.*sin(de1).*i7+d4.*sin(de1).*i8);
            ayx2=sin(de).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*sin(de1).*i5+d2.*sin(de1).*i6-d3.*cos(de1).*i7-d4.*cos(de1).*i8);
            ayy2=-cos(de).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*sin(de1).*i5+d2.*sin(de1).*i6-d3.*cos(de1).*i7-d4.*cos(de1).*i8);
            azx2=sin(de).*(c1.*Z4+c2.*Z12+c3.*Z14).*(d7.*i5+d8.*i6);
            azy2=-cos(de).*(c1.*Z4+c2.*Z12+c3.*Z14).*(d7.*i5+d8.*i6);
            % Next
            c1=-l.^2.*wce.^2./nper./ve;
            c2=-(n-l).*wce./ve;
            axx3=c1.*sin(de).*(Z2+c2.*Z1).*(d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12);
            axy3=-c1.*cos(de).*(Z2+c2.*Z1).*(d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12);
            ayx3=c1.*sin(de).*(Z2+c2.*Z1).*(d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12);
            ayy3=-c1.*cos(de).*(Z2+c2.*Z1).*(d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12);
            azx3=c1.*sin(de).*(Z12+c2.*Z4).*(d7.*i9+d8.*i10);
            azy3=-c1.*cos(de).*(Z12+c2.*Z4).*(d7.*i9+d8.*i10);
            % Next
            c1=1i.*l.*wce;
            axx4=c1.*cos(de).*(Z2+c2.*Z1).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
            axy4=c1.*sin(de).*(Z2+c2.*Z1).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
            axz4=1i.*nper.*(Z12-(n-l).*wce./ve.*Z4).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
            ayx4=c1.*cos(de).*(Z2+c2.*Z1).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
            ayy4=c1.*sin(de).*(Z2+c2.*Z1).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
            ayz4=1i.*nper.*(Z12-(n-l).*wce./ve.*Z4).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
            azx4=c1.*cos(de).*(Z12+c2.*Z4).*(d7.*i13+d8.*i14);
            azy4=c1.*sin(de).*(Z12+c2.*Z4).*(d7.*i13+d8.*i14);
            azz4=1i.*nper.*(Z13./ve-n0z.*Z18).*(d7.*i13+d8.*i14);
            % Sum
            sxx1=30./w0./nz./n1z.*(axx1+axx2+axx3+axx4);
            sxy1=30./w0./nz./n1z.*(axy1+axy2+axy3+axy4);
            sxz1=30./w0./nz./n1z.*ve.*(axz1+axz4);
            syx1=30./w0./nz./n1z.*(ayx1+ayx2+ayx3+ayx4);
            syy1=30./w0./nz./n1z.*(ayy1+ayy2+ayy3+ayy4);
            syz1=30./w0./nz./n1z.*ve.*(ayz1+ayz4);
            szx1=30./w0./nz./n1z.*(azx1+azx2+azx3+azx4);
            szy1=30./w0./nz./n1z.*(azy1+azy2+azy3+azy4);
            szz1=30./w0./nz./n1z.*ve.*(azz1+azz4);
    % Part 2
            d5=n.*wce./n1per./ve;
            c1=b2.*l.*wce;
            c2=-(n-l).*nper.^2.*wce;
            c3=1i.*l.*(1+(n-l).*wce).*wce.^2./ve.^2;
            axx1=((c1.*cos(de)+c2.*sin(de)).*Z4+c3.*cos(de).*Z5)./nper.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            axy1=((c1.*sin(de)-c2.*cos(de)).*Z4+c3.*sin(de).*Z5)./nper.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            axz1=(b2.*ve.*Z13+1i.*wce.*(1+(n-l).*wce)./ve.*zl.*Z5).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            ayx1=((c1.*cos(de)+c2.*sin(de)).*Z4+c3.*cos(de).*Z5)./nper.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            ayy1=((c1.*sin(de)-c2.*cos(de)).*Z4+c3.*sin(de).*Z5)./nper.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            ayz1=(b2.*ve.*Z13+1i.*wce.*(1+(n-l).*wce)./ve.*zl.*Z5).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            azx1=((c1.*cos(de)+c2.*sin(de)).*Z13+c3.*cos(de).*Z15)./nper.*i1;
            azy1=((c1.*sin(de)-c2.*cos(de)).*Z13+c3.*sin(de).*Z15)./nper.*i1;
            azz1=(b2.*ve.*Z18+1i.*wce.*(1+(n-l).*wce)./ve.*zl.*Z15).*i1;
            % Next
            c1=1i.*b2;
            c2=-(1+(n-l).*wce).*wce./ve.^2;
            axx2=ve.*sin(de).*(c1.*Z4+c2.*Z5).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
            axy2=-ve.*cos(de).*(c1.*Z4+c2.*Z5).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
            ayx2=ve.*sin(de).*(c1.*Z4+c2.*Z5).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
            ayy2=-ve.*cos(de).*(c1.*Z4+c2.*Z5).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
            azx2=ve.*sin(de).*(c1.*Z13+c2.*Z15).*i5;
            azy2=-ve.*cos(de).*(c1.*Z13+c2.*Z15).*i5;
            % Next
            c1=l.^2.*(n-l).*wce.^3./nper./ve.^2;
            axx3=c1.*sin(de).*Z4.*(d5.*cos(de1).*i9-1i.*sin(de1).*i11);
            axy3=-c1.*cos(de).*Z4.*(d5.*cos(de1).*i9-1i.*sin(de1).*i11);
            ayx3=c1.*sin(de).*Z4.*(d5.*sin(de1).*i9+1i.*cos(de1).*i11);
            ayy3=-c1.*cos(de).*Z4.*(d5.*sin(de1).*i9+1i.*cos(de1).*i11);
            azx3=c1.*sin(de).*Z13.*i9;
            azy3=-c1.*cos(de).*Z13.*i9;
            % Next
            c1=-1i.*l.*(n-l).*wce.^2./ve;
            axx4=c1.*cos(de).*Z4.*(d5.*cos(de1).*i13-1i.*sin(de1).*i15);
            axy4=c1.*sin(de).*Z4.*(d5.*cos(de1).*i13-1i.*sin(de1).*i15);
            axz4=-1i.*(n-l).*nper.*wce.*Z13.*(d5.*cos(de1).*i13-1i.*sin(de1).*i15);
            ayx4=c1.*cos(de).*Z4.*(d5.*sin(de1).*i13+1i.*cos(de1).*i15);
            ayy4=c1.*sin(de).*Z4.*(d5.*sin(de1).*i13+1i.*cos(de1).*i15);
            ayz4=-1i.*(n-l).*nper.*wce.*Z13.*(d5.*sin(de1).*i13+1i.*cos(de1).*i15);
            azx4=c1.*cos(de).*Z13.*i13;
            azy4=c1.*sin(de).*Z13.*i13;
            azz4=-1i.*(n-l).*nper.*wce.*Z18.*i13;
            % Sum
            sxx2=30./w0./nz./n1z.*E0zc.*(axx1+axx2+axx3+axx4);
            sxy2=30./w0./nz./n1z.*E0zc.*(axy1+axy2+axy3+axy4);
            sxz2=30./w0./nz./n1z.*E0zc.*(axz1+axz4);
            syx2=30./w0./nz./n1z.*E0zc.*(ayx1+ayx2+ayx3+ayx4);
            syy2=30./w0./nz./n1z.*E0zc.*(ayy1+ayy2+ayy3+ayy4);
            syz2=30./w0./nz./n1z.*E0zc.*(ayz1+ayz4);
            szx2=30./w0./nz./n1z.*E0zc.*(azx1+azx2+azx3+azx4);
            szy2=30./w0./nz./n1z.*E0zc.*(azy1+azy2+azy3+azy4);
            szz2=30./w0./nz./n1z.*E0zc.*(azz1+azz4);
    % Part 3
            d6=l.*wce./nper./ve;
            axx1=d6.*cos(de).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            axy1=d6.*sin(de).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            ayx1=d6.*cos(de).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            ayy1=d6.*sin(de).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            axx2=1i.*sin(de).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
            axy2=-1i.*cos(de).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
            ayx2=1i.*sin(de).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
            ayy2=-1i.*cos(de).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
            % Sum
            sxx3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z2.*(axx1+axx2);
            sxy3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z2.*(axy1+axy2);
            sxz3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z12.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            syx3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z2.*(ayx1+ayx2);
            syy3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z2.*(ayy1+ayy2);
            syz3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z12.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            szx3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z12.*(d6.*cos(de).*i1+1i.*sin(de).*i5);
            szy3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*ve.*Z12.*(d6.*sin(de).*i1-1i.*cos(de).*i5);
            szz3=30./w0./nz./n1z.*(nx.*E0yc-ny.*E0xc).*n0z.*ve.*(zl0.*Z13-Z18).*i1;
            % Next
            d1=-1i.*E0yc.*n.*(n-l).*wce.^2./n0x./n1per./ve.^2;
            d2=E0xc.*n.*wce./n1per./ve;
            d3=-E0yc.*(n-l).*wce./n0x./ve;
            d4=-1i.*E0xc;
            d7=-1i.*E0yc.*(n-l).*wce./n0x./ve;
            axx1=d6.*cos(de).*(d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12);
            axy1=d6.*sin(de).*(d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12);
            axz1=d1.*cos(de1).*i9+d2.*cos(de1).*i10+d3.*sin(de1).*i11+d4.*sin(de1).*i12;
            ayx1=d6.*cos(de).*(d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12);
            ayy1=d6.*sin(de).*(d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12);
            ayz1=d1.*sin(de1).*i9+d2.*sin(de1).*i10-d3.*cos(de1).*i11-d4.*cos(de1).*i12;
            azx1=d6.*cos(de).*(E0xc.*i10+d7.*i9);
            azy1=d6.*sin(de).*(E0xc.*i10+d7.*i9);
            azz1=E0xc.*i10+d7.*i9;
            axx2=1i.*sin(de).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
            axy2=-1i.*cos(de).*(d1.*cos(de1).*i13+d2.*cos(de1).*i14+d3.*sin(de1).*i15+d4.*sin(de1).*i16);
            ayx2=1i.*sin(de).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
            ayy2=-1i.*cos(de).*(d1.*sin(de1).*i13+d2.*sin(de1).*i14-d3.*cos(de1).*i15-d4.*cos(de1).*i16);
            azx2=1i.*sin(de).*(E0xc.*i14+d7.*i13);
            azy2=-1i.*cos(de).*(E0xc.*i14+d7.*i13);
            % Sum
            sxx4=30./w0./nz./n1z.*1i.*l.*wce.*(Z2-(n-l).*wce./ve.*Z1).*(axx1+axx2);
            sxy4=30./w0./nz./n1z.*1i.*l.*wce.*(Z2-(n-l).*wce./ve.*Z1).*(axy1+axy2);
            sxz4=30./w0./nz./n1z.*1i.*l.*wce.*(Z12-(n-l).*wce./ve.*Z4).*axz1;
            syx4=30./w0./nz./n1z.*1i.*l.*wce.*(Z2-(n-l).*wce./ve.*Z1).*(ayx1+ayx2);
            syy4=30./w0./nz./n1z.*1i.*l.*wce.*(Z2-(n-l).*wce./ve.*Z1).*(ayy1+ayy2);
            syz4=30./w0./nz./n1z.*1i.*l.*wce.*(Z12-(n-l).*wce./ve.*Z4).*ayz1;
            szx4=30./w0./nz./n1z.*1i.*l.*wce.*(Z12-(n-l).*wce./ve.*Z4).*(azx1+azx2);
            szy4=30./w0./nz./n1z.*1i.*l.*wce.*(Z12-(n-l).*wce./ve.*Z4).*(azy1+azy2);
            szz4=30./w0./nz./n1z.*1i.*l.*wce.*(Z13./ve-n0z.*Z18).*azz1;
            % Next
            axx1=d6.*cos(de).*(d5.*cos(de1).*i10-1i.*sin(de1).*i12);
            axy1=d6.*sin(de).*(d5.*cos(de1).*i10-1i.*sin(de1).*i12);
            ayx1=d6.*cos(de).*(d5.*sin(de1).*i10+1i.*cos(de1).*i12);
            ayy1=d6.*sin(de).*(d5.*sin(de1).*i10+1i.*cos(de1).*i12);
            axx2=1i.*sin(de).*(d5.*cos(de1).*i14-1i.*sin(de1).*i16);
            axy2=-1i.*cos(de).*(d5.*cos(de1).*i14-1i.*sin(de1).*i16);
            ayx2=1i.*sin(de).*(d5.*sin(de1).*i14+1i.*cos(de1).*i16);
            ayy2=-1i.*cos(de).*(d5.*sin(de1).*i14+1i.*cos(de1).*i16);
            % Sum
            sxx5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z4.*(axx1+axx2);
            sxy5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z4.*(axy1+axy2);
            sxz5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z13.*(d5.*cos(de1).*i10-1i.*sin(de1).*i12);
            syx5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z4.*(ayx1+ayx2);
            syy5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z4.*(ayy1+ayy2);
            syz5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z13.*(d5.*sin(de1).*i10+1i.*cos(de1).*i12);
            szx5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z13.*(d6.*cos(de).*i10+1i.*sin(de).*i14);
            szy5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z13.*(d6.*sin(de).*i10-1i.*cos(de).*i14);
            szz5=30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0zc.*Z18.*i10;
            % Next
            axx1=d6.*cos(de).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            axy1=d6.*sin(de).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            ayx1=d6.*cos(de).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            ayy1=d6.*sin(de).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            axx2=1i.*sin(de).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
            axy2=-1i.*cos(de).*(d5.*cos(de1).*i5-1i.*sin(de1).*i7);
            ayx2=1i.*sin(de).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
            ayy2=-1i.*cos(de).*(d5.*sin(de1).*i5+1i.*cos(de1).*i7);
            % Sum
            sxx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z1.*(axx1+axx2);
            sxy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z1.*(axy1+axy2);
            sxz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z4.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            syx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z1.*(ayx1+ayx2);
            syy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z1.*(ayy1+ayy2);
            syz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z4.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            szx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z4.*(d6.*cos(de).*i1+1i.*sin(de).*i5);
            szy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z4.*(d6.*sin(de).*i1-1i.*cos(de).*i5);
            szz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0yc.*Z13.*i1;
    % Part 4
            d1=n.*l.*wce.^2./nper./n1per./ve.^2;
            d2=1i.*n.*wce./n1per./ve;
            d3=-1i.*l.*wce./nper./ve;
            c1=b3.*E0zc-2i.*(n-l).*nz.*wce.^2./n0x./ve.^2.*E0xc;
            c2=(-(n-l).*b3./n0x.*E0xc-l.*n0x.*E0yc).*wce./ve;
            c3=-(-2i.*(n-l).*wce.^2./n0x./ve.^2.*E0xc-n0x.*E0yc);
            c4=-2i.*wce./ve.*E0zc;
            c5=2i.*nz.*wce./ve.*E0zc;
            c6=1i.*wce.*((1+(n-l).*wce).*E0zc-(n-l).*n0z.*wce./n0x.*E0xc).*nz./n0z./ve.^2;
            c7=c1.*Z9+c2.*Z6+c3.*Z7+c4.*Z10+c5.*Z11+c6.*Z8;
            c8=b4.*ve;
            c9=-l.*n0x.*wce.*E0yc;
            c10=1i.*wce.*(w-l.*wce).*((1+(n-l).*wce).*E0zc-(n-l).*n0z.*wce./n0x.*E0xc)./n0z./ve.^2;
            c11=c8.*(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9)+c9.*Z9+c10.*Z8;
            c12=(n0x.*ny+2i.*w.*wce./ve.^2).*E0zc-2i.*(n-l).*nz.*wce.^2./n0x./ve.^2.*E0xc;
            c13=(-b4.*(n-l)./n0x.*E0xc-l.*n0x.*E0yc).*wce./ve;
            c14=-1i.*wce.*((1+(n-l).*wce).*E0zc-(n-l).*n0z.*wce./n0x.*E0xc).*nz./n0z./ve.^2;
            axx1=c7.*(d1.*cos(de).*cos(de1).*i1+d2.*sin(de).*cos(de1).*i5+d3.*cos(de).*sin(de1).*i3+sin(de).*sin(de1).*i7);
            axy1=c7.*(d1.*sin(de).*cos(de1).*i1-d2.*cos(de).*cos(de1).*i5+d3.*sin(de).*sin(de1).*i3-cos(de).*sin(de1).*i7);
            axz1=c11.*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            ayx1=c7.*(d1.*cos(de).*sin(de1).*i1+d2.*sin(de).*sin(de1).*i5-d3.*cos(de).*cos(de1).*i3-sin(de).*cos(de1).*i7);
            ayy1=c7.*(d1.*sin(de).*sin(de1).*i1-d2.*cos(de).*sin(de1).*i5-d3.*sin(de).*cos(de1).*i3+cos(de).*cos(de1).*i7);
            ayz1=c11.*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            azx1=(c12.*Z11+c2.*Z9+c3.*Z10-c6.*Z16).*(d6.*cos(de).*i1+1i.*sin(de).*i5);
            azy1=(c12.*Z11+c2.*Z9+c3.*Z10-c6.*Z16).*(d6.*sin(de).*i1-1i.*cos(de).*i5);
            azz1=(b4.*E0zc.*Z17+c13.*Z11+c14.*zl.*Z16).*i1;
            % Next
            c1=1i.*b3;
            c2=2.*wce./ve;
            c3=-nz.*wce./ve;
            c4=E0yc.*(c1.*Z6+c2.*Z7+c3.*(2.*Z9+Z8));
            axx2=c4.*(d1.*cos(de).*cos(de1).*i2+d2.*sin(de).*cos(de1).*i6+d3.*cos(de).*sin(de1).*i4+sin(de).*sin(de1).*i8);
            axy2=c4.*(d1.*sin(de).*cos(de1).*i2-d2.*cos(de).*cos(de1).*i6+d3.*sin(de).*sin(de1).*i4-cos(de).*sin(de1).*i8);
            axz2=ve.*E0yc.*(1i.*b4.*Z9-wce.*(w-l.*wce)./ve.^2.*Z8).*(d5.*cos(de1).*i2-1i.*sin(de1).*i4);
            ayx2=c4.*(d1.*cos(de).*sin(de1).*i2+d2.*sin(de).*sin(de1).*i6-d3.*cos(de).*cos(de1).*i4-sin(de).*cos(de1).*i8);
            ayy2=c4.*(d1.*sin(de).*sin(de1).*i2-d2.*cos(de).*sin(de1).*i6-d3.*sin(de).*cos(de1).*i4+cos(de).*cos(de1).*i8);
            ayz2=ve.*E0yc.*(1i.*b4.*Z9-wce.*(w-l.*wce)./ve.^2.*Z8).*(d5.*sin(de1).*i2+1i.*cos(de1).*i4);
            azx2=E0yc.*(c1.*Z9+c2.*Z10+c3.*(2.*Z11-Z16)).*(d6.*cos(de).*i2+1i.*sin(de).*i6);
            azy2=E0yc.*(c1.*Z9+c2.*Z10+c3.*(2.*Z11-Z16)).*(d6.*sin(de).*i2-1i.*cos(de).*i6);
            azz2=E0yc.*(1i.*b4.*Z11+nz.*wce./ve.*zl.*Z16).*i2;
            % Next
            c1=-(n-l).^2.*wce.^2./n0x./ve.^2.*E0yc.*(Z7-l.*wce./ve.*Z6);
            axx3=c1.*(d1.*cos(de).*cos(de1).*i9+d2.*sin(de).*cos(de1).*i13+d3.*cos(de).*sin(de1).*i11+sin(de).*sin(de1).*i15);
            axy3=c1.*(d1.*sin(de).*cos(de1).*i9-d2.*cos(de).*cos(de1).*i13+d3.*sin(de).*sin(de1).*i11-cos(de).*sin(de1).*i15);
            axz3=E0yc.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z9.*(d5.*cos(de1).*i9-1i.*sin(de1).*i11);
            ayx3=c1.*(d1.*cos(de).*sin(de1).*i9+d2.*sin(de).*sin(de1).*i13-d3.*cos(de).*cos(de1).*i11-sin(de).*cos(de1).*i15);
            ayy3=c1.*(d1.*sin(de).*sin(de1).*i9-d2.*cos(de).*sin(de1).*i13-d3.*sin(de).*cos(de1).*i11+cos(de).*cos(de1).*i15);
            ayz3=E0yc.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z9.*(d5.*sin(de1).*i9+1i.*cos(de1).*i11);
            azx3=-(n-l).^2.*wce.^2./n0x./ve.^2.*E0yc.*(Z10-l.*wce./ve.*Z9).*(d6.*cos(de).*i9+1i.*sin(de).*i13);
            azy3=-(n-l).^2.*wce.^2./n0x./ve.^2.*E0yc.*(Z10-l.*wce./ve.*Z9).*(d6.*sin(de).*i9-1i.*cos(de).*i13);
            azz3=l.*(n-l).^2.*wce.^3./n0x./ve.^3.*E0yc.*Z11.*i9;
            % Next
            c1=1i.*n0x.*E0zc;
            c2=-1i.*l.*n0x.*wce./ve.*E0zc;
            c3=-1i.*(n-l).*wce./ve.*E0xc;
            c4=1i.*l.*(n-l).*wce.^2./ve.^2.*E0xc;
            c5=c1.*Z10+c2.*Z9+c3.*Z7+c4.*Z6;
            c6=-1i.*n0x.*w./ve.*E0zc;
            c7=1i.*nz.*n0x.*E0zc;
            axx4=c5.*(d1.*cos(de).*cos(de1).*i10+d2.*sin(de).*cos(de1).*i14+d3.*cos(de).*sin(de1).*i12+sin(de).*sin(de1).*i16);
            axy4=c5.*(d1.*sin(de).*cos(de1).*i10-d2.*cos(de).*cos(de1).*i14+d3.*sin(de).*sin(de1).*i12-cos(de).*sin(de1).*i16);
            axz4=-1i.*l.*n0x.*wce.*(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*(d5.*cos(de1).*i10-1i.*sin(de1).*i12);
            ayx4=c5.*(d1.*cos(de).*sin(de1).*i10+d2.*sin(de).*sin(de1).*i14-d3.*cos(de).*cos(de1).*i12-sin(de).*cos(de1).*i16);
            ayy4=c5.*(d1.*sin(de).*sin(de1).*i10-d2.*cos(de).*sin(de1).*i14-d3.*sin(de).*cos(de1).*i12+cos(de).*cos(de1).*i16);
            ayz4=-1i.*l.*n0x.*wce.*(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*(d5.*sin(de1).*i10+1i.*cos(de1).*i12);
            azx4=(c6.*Z11+c7.*Z17+c3.*Z10+c4.*Z9).*(d6.*cos(de).*i10+1i.*sin(de).*i14);
            azy4=(c6.*Z11+c7.*Z17+c3.*Z10+c4.*Z9).*(d6.*sin(de).*i10-1i.*cos(de).*i14);
            azz4=-1i.*l.*n0x.*wce./ve.*(E0zc.*Z17-E0xc.*(n-l).*wce./n0x./ve.*Z11).*i10;
            % Sum
            sxx7=30./w0./n0z./n1z.*ve./w.*(axx1+axx2+axx3+axx4);
            sxy7=30./w0./n0z./n1z.*ve./w.*(axy1+axy2+axy3+axy4);
            sxz7=30./w0./n0z./n1z./w.*(axz1+axz2+axz3+axz4);
            syx7=30./w0./n0z./n1z.*ve./w.*(ayx1+ayx2+ayx3+ayx4);
            syy7=30./w0./n0z./n1z.*ve./w.*(ayy1+ayy2+ayy3+ayy4);
            syz7=30./w0./n0z./n1z./w.*(ayz1+ayz2+ayz3+ayz4);
            szx7=30./w0./n0z./n1z.*ve./w.*(azx1+azx2+azx3+azx4);
            szy7=30./w0./n0z./n1z.*ve./w.*(azy1+azy2+azy3+azy4);
            szz7=30./w0./n0z./n1z.*ve./w.*(azz1+azz2+azz3+azz4);
    % Part 5
            axy1=(E0zc.*Z10-E0xc.*(n-l).*wce./n0x./ve.*Z7).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            ayy1=(E0zc.*Z10-E0xc.*(n-l).*wce./n0x./ve.*Z7).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            azy1=(E0zc.*((w-l.*wce)./ve.*Z11-nz.*Z17)-E0xc.*(n-l).*wce./n0x./ve.*(-Z10)).*i1;
            axy2=1i.*E0yc.*Z7.*(d5.*cos(de1).*i2-1i.*sin(de1).*i4);
            ayy2=1i.*E0yc.*Z7.*(d5.*sin(de1).*i2+1i.*cos(de1).*i4);
            azy2=-1i.*E0yc.*Z10.*i2;
            % Sum
            sxy8=-30./w0./n0z./n1z.*n0x.*ve./w.*(axy1+axy2);
            syy8=-30./w0./n0z./n1z.*n0x.*ve./w.*(ayy1+ayy2);
            szy8=30./w0./n0z./n1z.*n0x.*ve./w.*(azy1+azy2);
            % Next
            c1=E0zc.*(-Z10+l.*wce./ve.*Z9)-E0xc.*(n-l).*wce./n0x./ve.*(-Z7+l.*wce./ve.*Z6);
            axx1=c1.*(-d1.*sin(de).*cos(de1).*i9+d2.*cos(de).*cos(de1).*i13-d3.*sin(de).*sin(de1).*i11+cos(de).*sin(de1).*i15);
            axy1=c1.*(d1.*cos(de).*cos(de1).*i9+d2.*sin(de).*cos(de1).*i13+d3.*cos(de).*sin(de1).*i11+sin(de).*sin(de1).*i15);
            axz1=(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*(d5.*cos(de1).*i13-1i.*sin(de1).*i15);
            ayx1=c1.*(-d1.*sin(de).*sin(de1).*i9+d2.*cos(de).*sin(de1).*i13+d3.*sin(de).*cos(de1).*i11-cos(de).*cos(de1).*i15);
            ayy1=c1.*(d1.*cos(de).*sin(de1).*i9+d2.*sin(de).*sin(de1).*i13-d3.*cos(de).*cos(de1).*i11-sin(de).*cos(de1).*i15);
            ayz1=(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*(d5.*sin(de1).*i13+1i.*cos(de1).*i15);
            azx1=(E0zc.*(w./ve.*Z11-nz.*Z17)-E0xc.*(n-l).*wce./n0x./ve.*(-Z10+l.*wce./ve.*Z9)).*(d6.*sin(de).*i9-1i.*cos(de).*i13);
            azy1=(E0zc.*(w./ve.*Z11-nz.*Z17)-E0xc.*(n-l).*wce./n0x./ve.*(-Z10+l.*wce./ve.*Z9)).*(-d6.*cos(de).*i9-1i.*sin(de).*i13);
            azz1=(E0zc.*Z17-E0xc.*(n-l).*wce./n0x./ve.*Z11).*i13;
            c1=1i.*E0yc.*(-Z7+l.*wce./ve.*Z6);
            axx2=c1.*(-d1.*sin(de).*cos(de1).*i10+d2.*cos(de).*cos(de1).*i14-d3.*sin(de).*sin(de1).*i12+cos(de).*sin(de1).*i16);
            axy2=c1.*(d1.*cos(de).*cos(de1).*i10+d2.*sin(de).*cos(de1).*i14+d3.*cos(de).*sin(de1).*i12+sin(de).*sin(de1).*i16);
            axz2=1i.*E0yc.*Z9.*(d5.*cos(de1).*i14-1i.*sin(de1).*i16);
            ayx2=c1.*(-d1.*sin(de).*sin(de1).*i10+d2.*cos(de).*sin(de1).*i14+d3.*sin(de).*cos(de1).*i12-cos(de).*cos(de1).*i16);
            ayy2=c1.*(d1.*cos(de).*sin(de1).*i10+d2.*sin(de).*sin(de1).*i14-d3.*cos(de).*cos(de1).*i12-sin(de).*cos(de1).*i16);
            ayz2=1i.*E0yc.*Z9.*(d5.*sin(de1).*i14+1i.*cos(de1).*i16);
            azx2=1i.*E0yc.*(-Z10+l.*wce./ve.*Z9).*(d6.*sin(de).*i10-1i.*cos(de).*i14);
            azy2=1i.*E0yc.*(-Z10+l.*wce./ve.*Z9).*(-d6.*cos(de).*i10-1i.*sin(de).*i14);
            azz2=1i.*E0yc.*Z11.*i14;
            % Sum
            sxx9=30./w0./n0z./n1z.*(n-l).*wce./w.*(axx1+axx2);
            sxy9=30./w0./n0z./n1z.*(n-l).*wce./w.*(axy1+axy2);
            sxz9=30./w0./n0z./n1z.*1i.*(n-l).*nper.*wce./w.*(axz1+axz2);
            syx9=30./w0./n0z./n1z.*(n-l).*wce./w.*(ayx1+ayx2);
            syy9=30./w0./n0z./n1z.*(n-l).*wce./w.*(ayy1+ayy2);
            syz9=30./w0./n0z./n1z.*1i.*(n-l).*nper.*wce./w.*(ayz1+ayz2);
            szx9=-30./w0./n0z./n1z.*(n-l).*wce./w.*(azx1+azx2);
            szy9=-30./w0./n0z./n1z.*(n-l).*wce./w.*(azy1+azy2);
            szz9=30./w0./n0z./n1z.*1i.*(n-l).*nper.*wce./w.*(azz1+azz2);
            % Next
            ax1=(E0zc.*Z9-E0xc.*(n-l).*wce./n0x./ve.*Z6).*(d5.*cos(de1).*i1-1i.*sin(de1).*i3);
            ay1=(E0zc.*Z9-E0xc.*(n-l).*wce./n0x./ve.*Z6).*(d5.*sin(de1).*i1+1i.*cos(de1).*i3);
            az1=(E0zc.*Z11-E0xc.*(n-l).*wce./n0x./ve.*Z9).*i1;
            ax2=1i.*E0yc.*Z6.*(d5.*cos(de1).*i2-1i.*sin(de1).*i4);
            ay2=1i.*E0yc.*Z6.*(d5.*sin(de1).*i2+1i.*cos(de1).*i4);
            az2=1i.*E0yc.*Z9.*i2;
            % Sum
            sxx10=30./w0./n0z./n1z.*(n-l).*ny.*wce./w.*(ax1+ax2);
            sxy10=-30./w0./n0z./n1z.*(n-l).*nx.*wce./w.*(ax1+ax2);
            syx10=30./w0./n0z./n1z.*(n-l).*ny.*wce./w.*(ay1+ay2);
            syy10=-30./w0./n0z./n1z.*(n-l).*nx.*wce./w.*(ay1+ay2);
            szx10=30./w0./n0z./n1z.*(n-l).*ny.*wce./w.*(az1+az2);
            szy10=-30./w0./n0z./n1z.*(n-l).*nx.*wce./w.*(az1+az2);
            % % %
            sxx=sxx1+sxx2+sxx3+sxx4+sxx5+sxx6+sxx7+sxx9+sxx10;
            sxy=sxy1+sxy2+sxy3+sxy4+sxy5+sxy6+sxy7+sxy8+sxy9+sxy10;
            sxz=sxz1+sxz2+sxz3+sxz4+sxz5+sxz6+sxz7+sxz9;
            syx=syx1+syx2+syx3+syx4+syx5+syx6+syx7+syx9+syx10;
            syy=syy1+syy2+syy3+syy4+syy5+syy6+syy7+syy8+syy9+syy10;
            syz=syz1+syz2+syz3+syz4+syz5+syz6+syz7+syz9;
            szx=szx1+szx2+szx3+szx4+szx5+szx6+szx7+szx9+szx10;
            szy=szy1+szy2+szy3+szy4+szy5+szy6+szy7+szy8+szy9+szy10;
            szz=szz1+szz2+szz3+szz4+szz5+szz6+szz7+szz9;
            SigNLxx=SigNLxx+exp(1i.*(n.*de1-l.*de)).*sxx;
            SigNLxy=SigNLxy+exp(1i.*(n.*de1-l.*de)).*sxy;
            SigNLxz=SigNLxz+exp(1i.*(n.*de1-l.*de)).*sxz;
            SigNLyx=SigNLyx+exp(1i.*(n.*de1-l.*de)).*syx;
            SigNLyy=SigNLyy+exp(1i.*(n.*de1-l.*de)).*syy;
            SigNLyz=SigNLyz+exp(1i.*(n.*de1-l.*de)).*syz;
            SigNLzx=SigNLzx+exp(1i.*(n.*de1-l.*de)).*szx;
            SigNLzy=SigNLzy+exp(1i.*(n.*de1-l.*de)).*szy;
            SigNLzz=SigNLzz+exp(1i.*(n.*de1-l.*de)).*szz;
        end
    end
    %{
    y=[SigNLxx SigNLxy SigNLxz
        SigNLyx SigNLyy SigNLyz
        SigNLzx SigNLzy SigNLzz];
    %}
    y = 0.096./81.9.*wpe.^2./w1./wce./ve.^2.*[SigNLxx SigNLxy SigNLxz
        SigNLyx SigNLyy SigNLyz
        SigNLzx SigNLzy SigNLzz];
    num = num+1;
end
y = y;