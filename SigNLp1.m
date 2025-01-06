%低频 非线性
%完成于2015/6/16
function y = SigNLp1(w)
global wpe wce ve;
global de1 de w0 nz n0z n1z nper n1per n0x n1x n1y;
global E0x E0y E0z;
w1=w-1;
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
        zn=(w-n.*wce)./nz./ve;
        zl1=(w1-l.*wce)./n1z./ve;
        zn0=(1-(n-l).*wce)./n0z./ve;
        Z1=(Z(zn)-Z(zl1))./(zn-zl1);
        Z2=n1z.*Z(zn)-nz.*Z(zl1);
        Z3=-(1+zn.*Z(zn))./(zn-zl1)-(Z(zn)-Z(zl1))./2./(zn-zl1).^2;
        Z4=(zn.*Z(zn)-zl1.*Z(zl1))./(zn-zl1);
        Z5=-2.*(1+zl1.*Z(zl1))./(zn-zl1)-(Z(zn)-Z(zl1))./(zn-zl1).^2;
        Z6=(Z(zn)-Z(zn0))./(zn-zn0);
        Z7=n0z.*Z(zn)-nz.*Z(zn0);
        Z8=-2.*(1+zn0.*Z(zn0))./(zn-zn0)-(Z(zn)-Z(zn0))./(zn-zn0).^2;
        Z9=(zn.*Z(zn)-zn0.*Z(zn0))./(zn-zn0);
        Z10=n0z.*(1+zn.*Z(zn))-nz.*(1+zn0.*Z(zn0));
        Z11=1+(zn.^2.*Z(zn)-zn0.^2.*Z(zn0))./(zn-zn0);
        Z12=n1z.*(1+zn.*Z(zn))-nz.*(1+zl1.*Z(zl1));
        Z13=1+(zn.^2.*Z(zn)-zl1.^2.*Z(zl1))./(zn-zl1);
        Z14=zn.*Z3;
        Z15=-2.*zl1.*(1+zl1.*Z(zl1))./(zn-zl1)-zn.*(Z(zn)-Z(zl1))./(zn-zl1).^2;
        Z16=-2.*zn0.*(1+zn0.*Z(zn0))./(zn-zn0)-zn.*(Z(zn)-Z(zn0))./(zn-zn0).^2;
        Z17=(zn.^2.*(1+zn.*Z(zn))-zn0.^2.*(1+zn0.*Z(zn0)))./(zn-zn0);
        Z18=(zn.^2.*(1+zn.*Z(zn))-zl1.^2.*(1+zl1.*Z(zl1)))./(zn-zl1);
        %贝塞尔函数积分
        p=nper.*ve./wce;
        p0=n0x.*ve./wce;
        p1=n1per.*ve./wce;
        f1=@(x) x.*exp(-x.^2).*besselj(n,p.*x).*besselj(n-l,p0.*x).*besselj(l,p1.*x);
        i1=integral(f1,0,inf);
        f2=@(x) x.^2.*exp(-x.^2).*besselj(n,p.*x) ...
            .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)).*besselj(l,p1.*x);
        i2=integral(f2,0,inf);
        f3=@(x) x.^2.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
            .*besselj(n-l,p0.*x).*besselj(l,p1.*x);
        i3=integral(f3,0,inf);
        f4=@(x) x.^3.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
            .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)).*besselj(l,p1.*x);
        i4=integral(f4,0,inf);
        f5=@(x) x.^2.*exp(-x.^2).*besselj(n,p.*x).*besselj(n-l,p0.*x) ...
            .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
        i5=integral(f5,0,inf);
        f6=@(x) x.^3.*exp(-x.^2).*besselj(n,p.*x) ...
            .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)) ...
            .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
        i6=integral(f6,0,inf);
        f7=@(x) x.^3.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
            .*besselj(n-l,p0.*x).*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
        i7=integral(f7,0,inf);
        f8=@(x) x.^4.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
            .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)) ...
            .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
        i8=integral(f8,0,inf);
        f9=@(x) exp(-x.^2)./x.*besselj(n,p.*x).*besselj(n-l,p0.*x).*besselj(l,p1.*x);
        i9=integral(f9,1e-39,inf);    %此积分会出警告
        f10=@(x) exp(-x.^2).*besselj(n,p.*x) ...
            .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)).*besselj(l,p1.*x);
        i10=integral(f10,0,inf);
        f11=@(x) exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
            .*besselj(n-l,p0.*x).*besselj(l,p1.*x);
        i11=integral(f11,0,inf);
        f12=@(x) x.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
            .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)).*besselj(l,p1.*x);
        i12=integral(f12,0,inf);
        f13=@(x) exp(-x.^2).*besselj(n,p.*x).*besselj(n-l,p0.*x) ...
            .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
        i13=integral(f13,0,inf);
        f14=@(x) x.*exp(-x.^2).*besselj(n,p.*x) ...
            .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)) ...
            .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
        i14=integral(f14,0,inf);
        f15=@(x) x.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
            .*besselj(n-l,p0.*x).*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
        i15=integral(f15,0,inf);
        f16=@(x) x.^2.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
            .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)) ...
            .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
        i16=integral(f16,0,inf);
%Part 1
        d1=E0x.*n.*(n-l).*wce.^2./nper./n0x./ve.^2;
        d2=-1i.*E0y.*n.*wce./nper./ve;
        d3=-1i.*E0x.*(n-l).*wce./n0x./ve;
        d4=-E0y;
        d7=E0x.*(n-l).*wce./n0x./ve;
        d8=-1i.*E0y;
        b1=-(n0x.*n1y+2i.*(n-l).*wce.^2./ve.^2);
        b2=-(n0x.*n1y+2i.*wce./ve.^2);
        c1=b1.*l.*wce;
        c2=(n-l).*n1per.^2.*wce;
        c3=-2i.*l.*wce.^2./ve;
        c4=n1per.^2.*ve;
        c5=2i.*l.*n0z.*wce.^2./ve;
        axx1=((c1.*cos(de1)+c2.*sin(de1)).*Z1+(c3.*cos(de1)+c4.*sin(de1)).*Z2 ...
            +c5.*cos(de1).*Z3)./n1per.*(d1.*cos(de).*i1+d2.*cos(de).*i2+d3.*sin(de).*i3+d4.*sin(de).*i4);
        axy1=((c1.*sin(de1)-c2.*cos(de1)).*Z1+(c3.*sin(de1)-c4.*cos(de1)).*Z2 ...
            +c5.*sin(de1).*Z3)./n1per.*(d1.*cos(de).*i1+d2.*cos(de).*i2+d3.*sin(de).*i3+d4.*sin(de).*i4);
        axz1=(b1.*Z4-2i.*wce./ve.*Z12+2i.*n0z.*wce./ve.*(-Z13+zl1./2.*Z5)) ...
            .*(d1.*cos(de).*i1+d2.*cos(de).*i2+d3.*sin(de).*i3+d4.*sin(de).*i4);
        ayx1=((c1.*cos(de1)+c2.*sin(de1)).*Z1+(c3.*cos(de1)+c4.*sin(de1)).*Z2 ...
            +c5.*cos(de1).*Z3)./n1per.*(d1.*sin(de).*i1+d2.*sin(de).*i2-d3.*cos(de).*i3-d4.*cos(de).*i4);
        ayy1=((c1.*sin(de1)-c2.*cos(de1)).*Z1+(c3.*sin(de1)-c4.*cos(de1)).*Z2 ...
            +c5.*sin(de1).*Z3)./n1per.*(d1.*sin(de).*i1+d2.*sin(de).*i2-d3.*cos(de).*i3-d4.*cos(de).*i4);
        ayz1=(b1.*Z4-2i.*wce./ve.*Z12+2i.*n0z.*wce./ve.*(-Z13+zl1./2.*Z5)) ...
            .*(d1.*sin(de).*i1+d2.*sin(de).*i2-d3.*cos(de).*i3-d4.*cos(de).*i4);
        azx1=((c1.*cos(de1)+c2.*sin(de1)).*Z4+(c3.*cos(de1)+c4.*sin(de1)).*Z12 ...
            +c5.*cos(de1).*Z14)./n1per.*(d7.*i1+d8.*i2);
        azy1=((c1.*sin(de1)-c2.*cos(de1)).*Z4+(c3.*sin(de1)-c4.*cos(de1)).*Z12 ...
            +c5.*sin(de1).*Z14)./n1per.*(d7.*i1+d8.*i2);
        azz1=(b2.*Z13+1i.*n0z.*wce./ve.*zl1.*Z15).*(d7.*i1+d8.*i2);
        %Next
        c1=1i.*b1.*ve;
        c2=2.*wce;
        c3=-2.*n0z.*wce;
        axx2=sin(de1).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*cos(de).*i5+d2.*cos(de).*i6+d3.*sin(de).*i7+d4.*sin(de).*i8);
        axy2=-cos(de1).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*cos(de).*i5+d2.*cos(de).*i6+d3.*sin(de).*i7+d4.*sin(de).*i8);
        ayx2=sin(de1).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*sin(de).*i5+d2.*sin(de).*i6-d3.*cos(de).*i7-d4.*cos(de).*i8);
        ayy2=-cos(de1).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*sin(de).*i5+d2.*sin(de).*i6-d3.*cos(de).*i7-d4.*cos(de).*i8);
        azx2=sin(de1).*(c1.*Z4+c2.*Z12+c3.*Z14).*(d7.*i5+d8.*i6);
        azy2=-cos(de1).*(c1.*Z4+c2.*Z12+c3.*Z14).*(d7.*i5+d8.*i6);
        %Next
        c1=-l.^2.*wce.^2./n1per./ve;
        c2=(n-l).*wce./ve;
        axx3=c1.*sin(de1).*(Z2+c2.*Z1).*(d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12);
        axy3=-c1.*cos(de1).*(Z2+c2.*Z1).*(d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12);
        ayx3=c1.*sin(de1).*(Z2+c2.*Z1).*(d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12);
        ayy3=-c1.*cos(de1).*(Z2+c2.*Z1).*(d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12);
        azx3=c1.*sin(de1).*(Z12+c2.*Z4).*(d7.*i9+d8.*i10);
        azy3=-c1.*cos(de1).*(Z12+c2.*Z4).*(d7.*i9+d8.*i10);
        %Next
        c1=1i.*l.*wce;
        axx4=c1.*cos(de1).*(Z2+c2.*Z1).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);    %此项为零？
        axy4=c1.*sin(de1).*(Z2+c2.*Z1).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);
        axz4=1i.*n1per.*(Z12+(n-l).*wce./ve.*Z4).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);
        ayx4=c1.*cos(de1).*(Z2+c2.*Z1).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
        ayy4=c1.*sin(de1).*(Z2+c2.*Z1).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
        ayz4=1i.*n1per.*(Z12+(n-l).*wce./ve.*Z4).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
        azx4=c1.*cos(de1).*(Z12+c2.*Z4).*(d7.*i13+d8.*i14);
        azy4=c1.*sin(de1).*(Z12+c2.*Z4).*(d7.*i13+d8.*i14);
        azz4=1i.*n1per.*(Z13./ve-n0z.*Z18).*(d7.*i13+d8.*i14);
        %Sum
        sxx1=30./w0./nz./n1z.*(axx1+axx2+axx3+axx4);
        sxy1=30./w0./nz./n1z.*(axy1+axy2+axy3+axy4);
        sxz1=30./w0./nz./n1z.*ve.*(axz1+axz4);
        syx1=30./w0./nz./n1z.*(ayx1+ayx2+ayx3+ayx4);
        syy1=30./w0./nz./n1z.*(ayy1+ayy2+ayy3+ayy4);
        syz1=30./w0./nz./n1z.*ve.*(ayz1+ayz4);
        szx1=30./w0./nz./n1z.*(azx1+azx2+azx3+azx4);
        szy1=30./w0./nz./n1z.*(azy1+azy2+azy3+azy4);
        szz1=30./w0./nz./n1z.*ve.*(azz1+azz4);
%Part 2
        d5=n.*wce./nper./ve;
        d6=-1i;
        c1=b2.*l.*wce;
        c2=(n-l).*n1per.^2.*wce;
        c3=1i.*l.*(1-(n-l).*wce).*wce.^2./ve.^2;
        axx1=((c1.*cos(de1)+c2.*sin(de1)).*Z4+c3.*cos(de1).*Z5)./n1per.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        axy1=((c1.*sin(de1)-c2.*cos(de1)).*Z4+c3.*sin(de1).*Z5)./n1per.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        axz1=(b2.*ve.*Z13+1i.*(1-(n-l).*wce).*wce./ve.*zl1.*Z5).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        ayx1=((c1.*cos(de1)+c2.*sin(de1)).*Z4+c3.*cos(de1).*Z5)./n1per.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        ayy1=((c1.*sin(de1)-c2.*cos(de1)).*Z4+c3.*sin(de1).*Z5)./n1per.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        ayz1=(b2.*ve.*Z13+1i.*(1-(n-l).*wce).*wce./ve.*zl1.*Z5).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        azx1=((c1.*cos(de1)+c2.*sin(de1)).*Z13+c3.*cos(de1).*Z15)./n1per.*i1;
        azy1=((c1.*sin(de1)-c2.*cos(de1)).*Z13+c3.*sin(de1).*Z15)./n1per.*i1;
        azz1=(b2.*ve.*Z18+1i.*wce.*(1-(n-l).*wce)./ve.*zl1.*Z15).*i1;
        %Next
        axx2=ve.*sin(de1).*(1i.*b2.*Z4-(1-(n-l).*wce).*wce./ve.^2.*Z5).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
        axy2=-ve.*cos(de1).*(1i.*b2.*Z4-(1-(n-l).*wce).*wce./ve.^2.*Z5).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
        ayx2=ve.*sin(de1).*(1i.*b2.*Z4-(1-(n-l).*wce).*wce./ve.^2.*Z5).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
        ayy2=-ve.*cos(de1).*(1i.*b2.*Z4-(1-(n-l).*wce).*wce./ve.^2.*Z5).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
        azx2=ve.*sin(de1).*(1i.*b2.*Z13-(1-(n-l).*wce).*wce./ve.^2.*Z15).*i5;
        azy2=-ve.*cos(de1).*(1i.*b2.*Z13-(1-(n-l).*wce).*wce./ve.^2.*Z15).*i5;
        %Next
        axx3=-l.^2.*(n-l).*wce.^3./n1per./ve.^2.*sin(de1).*Z4.*(d5.*cos(de).*i9+d6.*sin(de).*i11);
        axy3=l.^2.*(n-l).*wce.^3./n1per./ve.^2.*cos(de1).*Z4.*(d5.*cos(de).*i9+d6.*sin(de).*i11);
        ayx3=-l.^2.*(n-l).*wce.^3./n1per./ve.^2.*sin(de1).*Z4.*(d5.*sin(de).*i9-d6.*cos(de).*i11);
        ayy3=l.^2.*(n-l).*wce.^3./n1per./ve.^2.*cos(de1).*Z4.*(d5.*sin(de).*i9-d6.*cos(de).*i11);
        azx3=-l.^2.*(n-l).*wce.^3./n1per./ve.^2.*sin(de1).*Z13.*i9;
        azy3=l.^2.*(n-l).*wce.^3./n1per./ve.^2.*cos(de1).*Z13.*i9;
        %Next
        axx4=1i.*l.*(n-l).*wce.^2./ve.*cos(de1).*Z4.*(d5.*cos(de).*i13+d6.*sin(de).*i15);    %此项为零?
        axy4=1i.*l.*(n-l).*wce.^2./ve.*sin(de1).*Z4.*(d5.*cos(de).*i13+d6.*sin(de).*i15);
        axz4=1i.*(n-l).*n1per.*wce.*Z13.*(d5.*cos(de).*i13+d6.*sin(de).*i15);
        ayx4=1i.*l.*(n-l).*wce.^2./ve.*cos(de1).*Z4.*(d5.*sin(de).*i13-d6.*cos(de).*i15);
        ayy4=1i.*l.*(n-l).*wce.^2./ve.*sin(de1).*Z4.*(d5.*sin(de).*i13-d6.*cos(de).*i15);
        ayz4=1i.*(n-l).*n1per.*wce.*Z13.*(d5.*sin(de).*i13-d6.*cos(de).*i15);
        azx4=1i.*l.*(n-l).*wce.^2./ve.*cos(de1).*Z13.*i13;
        azy4=1i.*l.*(n-l).*wce.^2./ve.*sin(de1).*Z13.*i13;
        azz4=1i.*(n-l).*n1per.*wce.*Z18.*i13;
        %Sum
        sxx2=30./w0./nz./n1z.*E0z.*(axx1+axx2+axx3+axx4);
        sxy2=30./w0./nz./n1z.*E0z.*(axy1+axy2+axy3+axy4);
        sxz2=30./w0./nz./n1z.*E0z.*(axz1+axz4);
        syx2=30./w0./nz./n1z.*E0z.*(ayx1+ayx2+ayx3+ayx4);
        syy2=30./w0./nz./n1z.*E0z.*(ayy1+ayy2+ayy3+ayy4);
        syz2=30./w0./nz./n1z.*E0z.*(ayz1+ayz4);
        szx2=30./w0./nz./n1z.*E0z.*(azx1+azx2+azx3+azx4);
        szy2=30./w0./nz./n1z.*E0z.*(azy1+azy2+azy3+azy4);
        szz2=30./w0./nz./n1z.*E0z.*(azz1+azz4);
%Part 3
        d9=l.*wce./n1per./ve;
        axx1=l.*wce./n1per.*cos(de1).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        axy1=l.*wce./n1per.*sin(de1).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        ayx1=l.*wce./n1per.*cos(de1).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        ayy1=l.*wce./n1per.*sin(de1).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        axx2=1i.*ve.*sin(de1).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
        axy2=-1i.*ve.*cos(de1).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
        ayx2=1i.*ve.*sin(de1).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
        ayy2=-1i.*ve.*cos(de1).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
        %Sum
        sxx3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*Z2.*(axx1+axx2);
        sxy3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*Z2.*(axy1+axy2);
        sxz3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*ve.*Z12.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        syx3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*Z2.*(ayx1+ayx2);
        syy3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*Z2.*(ayy1+ayy2);
        syz3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*ve.*Z12.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        szx3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*ve.*Z12.*(d9.*cos(de1).*i1+1i.*sin(de1).*i5);
        szy3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*ve.*Z12.*(d9.*sin(de1).*i1-1i.*cos(de1).*i5);
        szz3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*n0z.*ve.*(zn0.*Z13-Z18).*i1;
        %Next
        d1=-1i.*E0y.*n.*(n-l).*wce.^2./nper./n0x./ve.^2;
        d2=E0x.*n.*wce./nper./ve;
        d3=-E0y.*(n-l).*wce./n0x./ve;
        d4=-1i.*E0x;
        d7=-1i.*E0y.*(n-l).*wce./n0x./ve;
        axx1=l.*wce./n1per./ve.*cos(de1).*(d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12);
        axy1=l.*wce./n1per./ve.*sin(de1).*(d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12);
        axz1=d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12;
        ayx1=l.*wce./n1per./ve.*cos(de1).*(d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12);
        ayy1=l.*wce./n1per./ve.*sin(de1).*(d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12);
        ayz1=d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12;
        azx1=d9.*cos(de1).*(d7.*i9+E0x.*i10);
        azy1=d9.*sin(de1).*(d7.*i9+E0x.*i10);
        azz1=d7.*i9+E0x.*i10;
        axx2=1i.*sin(de1).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);
        axy2=-1i.*cos(de1).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);
        ayx2=1i.*sin(de1).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
        ayy2=-1i.*cos(de1).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
        azx2=1i.*sin(de1).*(d7.*i13+E0x.*i14);
        azy2=-1i.*cos(de1).*(d7.*i13+E0x.*i14);
        %Sum
        sxx4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z2+(n-l).*wce./ve.*Z1).*(axx1+axx2);
        sxy4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z2+(n-l).*wce./ve.*Z1).*(axy1+axy2);
        sxz4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z12+(n-l).*wce./ve.*Z4).*axz1;
        syx4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z2+(n-l).*wce./ve.*Z1).*(ayx1+ayx2);
        syy4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z2+(n-l).*wce./ve.*Z1).*(ayy1+ayy2);
        syz4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z12+(n-l).*wce./ve.*Z4).*ayz1;
        szx4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z12+(n-l).*wce./ve.*Z4).*(azx1+azx2);
        szy4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z12+(n-l).*wce./ve.*Z4).*(azy1+azy2);
        szz4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z13./ve-n0z.*Z18).*azz1;
        %Next
        axx1=l.*wce./n1per./ve.*cos(de1).*(d5.*cos(de).*i10+d6.*sin(de).*i12);
        axy1=l.*wce./n1per./ve.*sin(de1).*(d5.*cos(de).*i10+d6.*sin(de).*i12);
        ayx1=l.*wce./n1per./ve.*cos(de1).*(d5.*sin(de).*i10-d6.*cos(de).*i12);
        ayy1=l.*wce./n1per./ve.*sin(de1).*(d5.*sin(de).*i10-d6.*cos(de).*i12);
        axx2=1i.*sin(de1).*(d5.*cos(de).*i14+d6.*sin(de).*i16);
        axy2=-1i.*cos(de1).*(d5.*cos(de).*i14+d6.*sin(de).*i16);
        ayx2=1i.*sin(de1).*(d5.*sin(de).*i14-d6.*cos(de).*i16);
        ayy2=-1i.*cos(de1).*(d5.*sin(de).*i14-d6.*cos(de).*i16);
        %Sum
        sxx5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z4.*(axx1+axx2);
        sxy5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z4.*(axy1+axy2);
        sxz5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z13.*(d5.*cos(de).*i10+d6.*sin(de).*i12);
        syx5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z4.*(ayx1+ayx2);
        syy5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z4.*(ayy1+ayy2);
        syz5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z13.*(d5.*sin(de).*i10-d6.*cos(de).*i12);
        szx5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z13.*(d9.*cos(de1).*i10+1i.*sin(de1).*i14);
        szy5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z13.*(d9.*sin(de1).*i10-1i.*cos(de1).*i14);
        szz5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z18.*i10;
        %Next
        axx1=l.*wce./n1per./ve.*cos(de1).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        axy1=l.*wce./n1per./ve.*sin(de1).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        ayx1=l.*wce./n1per./ve.*cos(de1).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        ayy1=l.*wce./n1per./ve.*sin(de1).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        axx2=1i.*sin(de1).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
        axy2=-1i.*cos(de1).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
        ayx2=1i.*sin(de1).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
        ayy2=-1i.*cos(de1).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
        %Sum
        sxx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z1.*(axx1+axx2);
        sxy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z1.*(axy1+axy2);
        sxz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z4.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        syx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z1.*(ayx1+ayx2);
        syy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z1.*(ayy1+ayy2);
        syz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z4.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        szx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z4.*(d9.*cos(de1).*i1+1i.*sin(de1).*i5);
        szy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z4.*(d9.*sin(de1).*i1-1i.*cos(de1).*i5);
        szz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z13.*i1;
%Part 4
        d1=n.*l.*wce.^2./nper./n1per./ve.^2;
        d2=1i.*n.*wce./nper./ve;
        d3=-1i.*l.*wce./n1per./ve;
        d4=1;
        b3=-2i.*(n-l).*wce.^2./n0x./ve.^2.*E0x-n0x.*E0y;
        c1=(n0x.*n1y-2i.*l.*wce.^2./ve.^2).*E0z-2i.*(n-l).*n1z.*wce.^2./n0x./ve.^2.*E0x;
        c2=((n-l).*n1y.*E0x+l.*b3).*wce./ve;
        c3=b3;
        c4=-2i.*wce./ve.*E0z;
        c5=-2i.*n1z.*wce./ve.*E0z;
        c6=1i.*wce.*((1-(n-l).*wce).*E0z+(n-l).*n0z./n0x.*wce.*E0x).*n1z./n0z./ve.^2;
        c7=c1.*Z9+c2.*Z6+c3.*Z7+c4.*Z10+c5.*Z11+c6.*Z8;
        c8=(n0x.*n1y-2i.*wce.*w1./ve.^2).*ve;
        c9=-l.*n0x.*wce.*E0y;
        c10=1i.*wce.*(w1-l.*wce).*((1-(n-l).*wce).*E0z+(n-l).*n0z.*wce./n0x.*E0x)./n0z./ve.^2;
        c11=c8.*(E0z.*Z11+E0x.*(n-l).*wce./n0x./ve.*Z9)+c9.*Z9+c10.*Z8;
        c12=(n0x.*n1y-2i.*w1.*wce./ve.^2).*E0z-2i.*(n-l).*n1z.*wce.^2./n0x./ve.^2.*E0x;
        axx1=c7.*(d1.*cos(de).*cos(de1).*i1+d2.*cos(de).*sin(de1).*i5+d3.*sin(de).*cos(de1).*i3+d4.*sin(de).*sin(de1).*i7);
        axy1=c7.*(d1.*cos(de).*sin(de1).*i1-d2.*cos(de).*cos(de1).*i5+d3.*sin(de).*sin(de1).*i3-d4.*sin(de).*cos(de1).*i7);
        axz1=c11.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        ayx1=c7.*(d1.*sin(de).*cos(de1).*i1+d2.*sin(de).*sin(de1).*i5-d3.*cos(de).*cos(de1).*i3-d4.*cos(de).*sin(de1).*i7);
        ayy1=c7.*(d1.*sin(de).*sin(de1).*i1-d2.*sin(de).*cos(de1).*i5-d3.*cos(de).*sin(de1).*i3+d4.*cos(de).*cos(de1).*i7);
        ayz1=c11.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        azx1=(c12.*Z11+c2.*Z9+c3.*Z10+c6.*Z16).*(d9.*cos(de1).*i1+1i.*sin(de1).*i5);
        azy1=(c12.*Z11+c2.*Z9+c3.*Z10+c6.*Z16).*(d9.*sin(de1).*i1-1i.*cos(de1).*i5);
        azz1=(c8.*(E0z.*Z17+E0x.*(n-l).*wce./n0x./ve.*Z11)+c9.*Z11+c10.*Z16).*i1;
        %Next
        c1=(-1i.*n0x.*n1y-2.*l.*wce.^2./ve.^2).*E0y;
        c2=-2.*wce./ve.*E0y;
        c3=n1z.*wce./ve.*E0y;
        c4=c1.*Z6+c2.*Z7+c3.*(-2.*Z9+Z8);
        axx2=c4.*(d1.*cos(de).*cos(de1).*i2+d2.*cos(de).*sin(de1).*i6+d3.*sin(de).*cos(de1).*i4+d4.*sin(de).*sin(de1).*i8);
        axy2=c4.*(d1.*cos(de).*sin(de1).*i2-d2.*cos(de).*cos(de1).*i6+d3.*sin(de).*sin(de1).*i4-d4.*sin(de).*cos(de1).*i8);
        axz2=ve.*E0y.*(-1i.*(n0x.*n1y-2i.*wce.*w1./ve.^2).*Z9+wce.*(w1-l.*wce)./ve.^2.*Z8).*(d5.*cos(de).*i2+d6.*sin(de).*i4);
        ayx2=c4.*(d1.*sin(de).*cos(de1).*i2+d2.*sin(de).*sin(de1).*i6-d3.*cos(de).*cos(de1).*i4-d4.*cos(de).*sin(de1).*i8);
        ayy2=c4.*(d1.*sin(de).*sin(de1).*i2-d2.*sin(de).*cos(de1).*i6-d3.*cos(de).*sin(de1).*i4+d4.*cos(de).*cos(de1).*i8);
        ayz2=ve.*E0y.*(-1i.*(n0x.*n1y-2i.*wce.*w1./ve.^2).*Z9+wce.*(w1-l.*wce)./ve.^2.*Z8).*(d5.*sin(de).*i2-d6.*cos(de).*i4);
        azx2=(c1.*Z9+c2.*Z10+c3.*(-2.*Z11+Z16)).*(d9.*cos(de1).*i2+1i.*sin(de1).*i6);
        azy2=(c1.*Z9+c2.*Z10+c3.*(-2.*Z11+Z16)).*(d9.*sin(de1).*i2-1i.*cos(de1).*i6);
        azz2=ve.*E0y.*(-1i.*(n0x.*n1y-2i.*wce.*w1./ve.^2).*Z11+wce.*(w1-l.*wce)./ve.^2.*Z16).*i2;
        %Next
        c1=(n-l).^2.*wce.^2./n0x./ve.^2.*E0y.*(Z7+l.*wce./ve.*Z6);
        axx3=c1.*(d1.*cos(de).*cos(de1).*i9+d2.*cos(de).*sin(de1).*i13+d3.*sin(de).*cos(de1).*i11+d4.*sin(de).*sin(de1).*i15);
        axy3=c1.*(d1.*cos(de).*sin(de1).*i9-d2.*cos(de).*cos(de1).*i13+d3.*sin(de).*sin(de1).*i11-d4.*sin(de).*cos(de1).*i15);
        axz3=E0y.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z9.*(d5.*cos(de).*i9+d6.*sin(de).*i11);
        ayx3=c1.*(d1.*sin(de).*cos(de1).*i9+d2.*sin(de).*sin(de1).*i13-d3.*cos(de).*cos(de1).*i11-d4.*cos(de).*sin(de1).*i15);
        ayy3=c1.*(d1.*sin(de).*sin(de1).*i9-d2.*sin(de).*cos(de1).*i13-d3.*cos(de).*sin(de1).*i11+d4.*cos(de).*cos(de1).*i15);
        ayz3=E0y.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z9.*(d5.*sin(de).*i9-d6.*cos(de).*i11);
        azx3=(n-l).^2.*wce.^2./n0x./ve.^2.*E0y.*(Z10+l.*wce./ve.*Z9).*(d9.*cos(de1).*i9+1i.*sin(de1).*i13);
        azy3=(n-l).^2.*wce.^2./n0x./ve.^2.*E0y.*(Z10+l.*wce./ve.*Z9).*(d9.*sin(de1).*i9-1i.*cos(de1).*i13);
        azz3=E0y.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z11.*i9;
        %Next
        c1=1i.*n0x.*E0z;
        c2=1i.*l.*n0x.*wce./ve.*E0z;
        c3=1i.*(n-l).*wce./ve.*E0x;
        c4=1i.*l.*(n-l).*wce.^2./ve.^2.*E0x;
        c5=c1.*Z10+c2.*Z9+c3.*Z7+c4.*Z6;
        c6=1i.*n0x.*w1./ve.*E0z;
        c7=-1i.*n0x.*n1z.*E0z;
        axx4=c5.*(d1.*cos(de).*cos(de1).*i10+d2.*cos(de).*sin(de1).*i14+d3.*sin(de).*cos(de1).*i12+d4.*sin(de).*sin(de1).*i16);
        axy4=c5.*(d1.*cos(de).*sin(de1).*i10-d2.*cos(de).*cos(de1).*i14+d3.*sin(de).*sin(de1).*i12-d4.*sin(de).*cos(de1).*i16);
        axz4=1i.*l.*wce.*(n0x.*E0z.*Z11+(n-l).*wce./ve.*E0x.*Z9).*(d5.*cos(de).*i10+d6.*sin(de).*i12);
        ayx4=c5.*(d1.*sin(de).*cos(de1).*i10+d2.*sin(de).*sin(de1).*i14-d3.*cos(de).*cos(de1).*i12-d4.*cos(de).*sin(de1).*i16);
        ayy4=c5.*(d1.*sin(de).*sin(de1).*i10-d2.*sin(de).*cos(de1).*i14-d3.*cos(de).*sin(de1).*i12+d4.*cos(de).*cos(de1).*i16);
        ayz4=1i.*l.*wce.*(n0x.*E0z.*Z11+(n-l).*wce./ve.*E0x.*Z9).*(d5.*sin(de).*i10-d6.*cos(de).*i12);
        azx4=(c6.*Z11+c7.*Z17+c3.*Z10+c4.*Z9).*(d9.*cos(de1).*i10+1i.*sin(de1).*i14);
        azy4=(c6.*Z11+c7.*Z17+c3.*Z10+c4.*Z9).*(d9.*sin(de1).*i10-1i.*cos(de1).*i14);
        azz4=1i.*l.*n0x.*wce.*(E0z.*Z17+E0x.*(n-l).*wce./n0x./ve.*Z11).*i10;
        %Sum
        sxx7=30./w0./nz./n0z.*ve./w1.*(axx1+axx2+axx3+axx4);
        sxy7=30./w0./nz./n0z.*ve./w1.*(axy1+axy2+axy3+axy4);
        sxz7=30./w0./nz./n0z./w1.*(axz1+axz2+axz3+axz4);
        syx7=30./w0./nz./n0z.*ve./w1.*(ayx1+ayx2+ayx3+ayx4);
        syy7=30./w0./nz./n0z.*ve./w1.*(ayy1+ayy2+ayy3+ayy4);
        syz7=30./w0./nz./n0z./w1.*(ayz1+ayz2+ayz3+ayz4);
        szx7=30./w0./nz./n0z.*ve./w1.*(azx1+azx2+azx3+azx4);
        szy7=30./w0./nz./n0z.*ve./w1.*(azy1+azy2+azy3+azy4);
        szz7=30./w0./nz./n0z./w1.*(azz1+azz2+azz3+azz4);
%Part 5
        axy1=(E0z.*Z10+E0x.*(n-l).*wce./n0x./ve.*Z7).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        ayy1=(E0z.*Z10+E0x.*(n-l).*wce./n0x./ve.*Z7).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        azy1=(E0z.*((w1-l.*wce)./ve.*Z11-n1z.*Z17)+E0x.*(n-l).*wce./n0x./ve.*Z10).*i1;
        axy2=-1i.*E0y.*Z7.*(d5.*cos(de).*i2+d6.*sin(de).*i4);
        ayy2=-1i.*E0y.*Z7.*(d5.*sin(de).*i2-d6.*cos(de).*i4);
        azy2=-1i.*E0y.*Z10.*i2;
        sxy8=30./w0./nz./n0z.*n0x.*ve./w1.*(axy1+axy2);
        syy8=30./w0./nz./n0z.*n0x.*ve./w1.*(ayy1+ayy2);
        szy8=30./w0./nz./n0z.*n0x.*ve./w1.*(azy1+azy2);
        %Next
        c1=E0z.*(Z10+l.*wce./ve.*Z9)+E0x.*(n-l).*wce./n0x./ve.*(Z7+l.*wce./ve.*Z6);
        axx1=c1.*(-d1.*cos(de).*sin(de1).*i9+d2.*cos(de).*cos(de1).*i13-d3.*sin(de).*sin(de1).*i11+d4.*sin(de).*cos(de1).*i15);
        axy1=c1.*(d1.*cos(de).*cos(de1).*i9+d2.*cos(de).*sin(de1).*i13+d3.*sin(de).*cos(de1).*i11+d4.*sin(de).*sin(de1).*i15);
        axz1=(E0z.*Z11+E0x.*(n-l).*wce./n0x./ve.*Z9).*(d5.*cos(de).*i13+d6.*sin(de).*i15);
        ayx1=c1.*(-d1.*sin(de).*sin(de1).*i9+d2.*sin(de).*cos(de1).*i13+d3.*cos(de).*sin(de1).*i11-d4.*cos(de).*cos(de1).*i15);
        ayy1=c1.*(d1.*sin(de).*cos(de1).*i9+d2.*sin(de).*sin(de1).*i13-d3.*cos(de).*cos(de1).*i11-d4.*cos(de).*sin(de1).*i15);
        ayz1=(E0z.*Z11+E0x.*(n-l).*wce./n0x./ve.*Z9).*(d5.*sin(de).*i13-d6.*cos(de).*i15);
        azx1=(E0z.*(w1./ve.*Z11-n1z.*Z17)+E0x.*(n-l).*wce./n0x./ve.*(Z10+l.*wce./ve.*Z9)).*(-d9.*sin(de1).*i9+1i.*cos(de1).*i13);
        azy1=(E0z.*(w1./ve.*Z11-n1z.*Z17)+E0x.*(n-l).*wce./n0x./ve.*(Z10+l.*wce./ve.*Z9)).*(d9.*cos(de1).*i9+1i.*sin(de1).*i13);
        azz1=(E0z.*Z17+E0x.*(n-l).*wce./n0x./ve.*Z11).*i13;
        c1=-1i.*E0y.*(Z7+l.*wce./ve.*Z6);
        axx2=c1.*(-d1.*cos(de).*sin(de1).*i10+d2.*cos(de).*cos(de1).*i14-d3.*sin(de).*sin(de1).*i12+d4.*sin(de).*cos(de1).*i16);
        axy2=c1.*(d1.*cos(de).*cos(de1).*i10+d2.*cos(de).*sin(de1).*i14+d3.*sin(de).*cos(de1).*i12+d4.*sin(de).*sin(de1).*i16);
        axz2=-1i.*E0y.*Z9.*(d5.*cos(de).*i14+d6.*sin(de).*i16);
        ayx2=c1.*(-d1.*sin(de).*sin(de1).*i10+d2.*sin(de).*cos(de1).*i14+d3.*cos(de).*sin(de1).*i12-d4.*cos(de).*cos(de1).*i16);
        ayy2=c1.*(d1.*sin(de).*cos(de1).*i10+d2.*sin(de).*sin(de1).*i14-d3.*cos(de).*cos(de1).*i12-d4.*cos(de).*sin(de1).*i16);
        ayz2=-1i.*E0y.*Z9.*(d5.*sin(de).*i14-d6.*cos(de).*i16);
        azx2=-1i.*E0y.*(Z10+l.*wce./ve.*Z9).*(-d9.*sin(de1).*i10+1i.*cos(de1).*i14);
        azy2=-1i.*E0y.*(Z10+l.*wce./ve.*Z9).*(d9.*cos(de1).*i10+1i.*sin(de1).*i14);
        azz2=-1i.*E0y.*Z11.*i14;
        %Sum
        sxx9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(axx1+axx2);
        sxy9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(axy1+axy2);
        sxz9=-30./w0./nz./n0z.*1i.*(n-l).*n1per.*wce./w1.*(axz1+axz2);
        syx9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(ayx1+ayx2);
        syy9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(ayy1+ayy2);
        syz9=-30./w0./nz./n0z.*1i.*(n-l).*n1per.*wce./w1.*(ayz1+ayz2);
        szx9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(azx1+azx2);
        szy9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(azy1+azy2);
        szz9=-30./w0./nz./n0z.*1i.*(n-l).*n1per.*wce./w1.*(azz1+azz2);
        %Next
        a1=(E0z.*Z9+E0x.*(n-l).*wce./n0x./ve.*Z6).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
        ay1=(E0z.*Z9+E0x.*(n-l).*wce./n0x./ve.*Z6).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
        az1=(E0z.*Z11+E0x.*(n-l).*wce./n0x./ve.*Z9).*i1;
        a2=-1i.*E0y.*Z6.*(d5.*cos(de).*i2+d6.*sin(de).*i4);
        ay2=-1i.*E0y.*Z6.*(d5.*sin(de).*i2-d6.*cos(de).*i4);
        az2=-1i.*E0y.*Z9.*i2;
        %Sum
        sxx10=-30./w0./nz./n0z.*(n-l).*n1y.*wce./w1.*(a1+a2);
        sxy10=30./w0./nz./n0z.*(n-l).*n1x.*wce./w1.*(a1+a2);
        syx10=-30./w0./nz./n0z.*(n-l).*n1y.*wce./w1.*(ay1+ay2);
        syy10=30./w0./nz./n0z.*(n-l).*n1x.*wce./w1.*(ay1+ay2);
        szx10=-30./w0./nz./n0z.*(n-l).*n1y.*wce./w1.*(az1+az2);
        szy10=30./w0./nz./n0z.*(n-l).*n1x.*wce./w1.*(az1+az2);
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
        SigNLxx=SigNLxx+exp(1i.*(n.*de-l.*de1)).*sxx;
        SigNLxy=SigNLxy+exp(1i.*(n.*de-l.*de1)).*sxy;
        SigNLxz=SigNLxz+exp(1i.*(n.*de-l.*de1)).*sxz;
        SigNLyx=SigNLyx+exp(1i.*(n.*de-l.*de1)).*syx;
        SigNLyy=SigNLyy+exp(1i.*(n.*de-l.*de1)).*syy;
        SigNLyz=SigNLyz+exp(1i.*(n.*de-l.*de1)).*syz;
        SigNLzx=SigNLzx+exp(1i.*(n.*de-l.*de1)).*szx;
        SigNLzy=SigNLzy+exp(1i.*(n.*de-l.*de1)).*szy;
        SigNLzz=SigNLzz+exp(1i.*(n.*de-l.*de1)).*szz;
    end
end
%{
y=[SigNLxx SigNLxy SigNLxz
    SigNLyx SigNLyy SigNLyz
    SigNLzx SigNLzy SigNLzz];
%}
y = 0.096./81.9.*wpe.^2./w./wce./ve.^2.*[SigNLxx SigNLxy SigNLxz
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
            zn=(w-n.*wce)./nz./ve;
            zl1=(w1-l.*wce)./n1z./ve;
            zn0=(1-(n-l).*wce)./n0z./ve;
            Z1=(Z(zn)-Z(zl1))./(zn-zl1);
            Z2=n1z.*Z(zn)-nz.*Z(zl1);
            Z3=-(1+zn.*Z(zn))./(zn-zl1)-(Z(zn)-Z(zl1))./2./(zn-zl1).^2;
            Z4=(zn.*Z(zn)-zl1.*Z(zl1))./(zn-zl1);
            Z5=-2.*(1+zl1.*Z(zl1))./(zn-zl1)-(Z(zn)-Z(zl1))./(zn-zl1).^2;
            Z6=(Z(zn)-Z(zn0))./(zn-zn0);
            Z7=n0z.*Z(zn)-nz.*Z(zn0);
            Z8=-2.*(1+zn0.*Z(zn0))./(zn-zn0)-(Z(zn)-Z(zn0))./(zn-zn0).^2;
            Z9=(zn.*Z(zn)-zn0.*Z(zn0))./(zn-zn0);
            Z10=n0z.*(1+zn.*Z(zn))-nz.*(1+zn0.*Z(zn0));
            Z11=1+(zn.^2.*Z(zn)-zn0.^2.*Z(zn0))./(zn-zn0);
            Z12=n1z.*(1+zn.*Z(zn))-nz.*(1+zl1.*Z(zl1));
            Z13=1+(zn.^2.*Z(zn)-zl1.^2.*Z(zl1))./(zn-zl1);
            Z14=zn.*Z3;
            Z15=-2.*zl1.*(1+zl1.*Z(zl1))./(zn-zl1)-zn.*(Z(zn)-Z(zl1))./(zn-zl1).^2;
            Z16=-2.*zn0.*(1+zn0.*Z(zn0))./(zn-zn0)-zn.*(Z(zn)-Z(zn0))./(zn-zn0).^2;
            Z17=(zn.^2.*(1+zn.*Z(zn))-zn0.^2.*(1+zn0.*Z(zn0)))./(zn-zn0);
            Z18=(zn.^2.*(1+zn.*Z(zn))-zl1.^2.*(1+zl1.*Z(zl1)))./(zn-zl1);
            %贝塞尔函数积分
            p=nper.*ve./wce;
            p0=n0x.*ve./wce;
            p1=n1per.*ve./wce;
            f1=@(x) x.*exp(-x.^2).*besselj(n,p.*x).*besselj(n-l,p0.*x).*besselj(l,p1.*x);
            i1=integral(f1,0,inf);
            f2=@(x) x.^2.*exp(-x.^2).*besselj(n,p.*x) ...
                .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)).*besselj(l,p1.*x);
            i2=integral(f2,0,inf);
            f3=@(x) x.^2.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
                .*besselj(n-l,p0.*x).*besselj(l,p1.*x);
            i3=integral(f3,0,inf);
            f4=@(x) x.^3.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
                .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)).*besselj(l,p1.*x);
            i4=integral(f4,0,inf);
            f5=@(x) x.^2.*exp(-x.^2).*besselj(n,p.*x).*besselj(n-l,p0.*x) ...
                .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
            i5=integral(f5,0,inf);
            f6=@(x) x.^3.*exp(-x.^2).*besselj(n,p.*x) ...
                .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)) ...
                .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
            i6=integral(f6,0,inf);
            f7=@(x) x.^3.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
                .*besselj(n-l,p0.*x).*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
            i7=integral(f7,0,inf);
            f8=@(x) x.^4.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
                .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)) ...
                .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
            i8=integral(f8,0,inf);
            f9=@(x) exp(-x.^2)./x.*besselj(n,p.*x).*besselj(n-l,p0.*x).*besselj(l,p1.*x);
            i9=integral(f9,1e-39,inf);    %此积分会出警告
            f10=@(x) exp(-x.^2).*besselj(n,p.*x) ...
                .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)).*besselj(l,p1.*x);
            i10=integral(f10,0,inf);
            f11=@(x) exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
                .*besselj(n-l,p0.*x).*besselj(l,p1.*x);
            i11=integral(f11,0,inf);
            f12=@(x) x.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
                .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)).*besselj(l,p1.*x);
            i12=integral(f12,0,inf);
            f13=@(x) exp(-x.^2).*besselj(n,p.*x).*besselj(n-l,p0.*x) ...
                .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
            i13=integral(f13,0,inf);
            f14=@(x) x.*exp(-x.^2).*besselj(n,p.*x) ...
                .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)) ...
                .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
            i14=integral(f14,0,inf);
            f15=@(x) x.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
                .*besselj(n-l,p0.*x).*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
            i15=integral(f15,0,inf);
            f16=@(x) x.^2.*exp(-x.^2).*(n.*besselj(n,p.*x)./(p.*x)-besselj(n+1,p.*x)) ...
                .*((n-l).*besselj(n-l,p0.*x)./(p0.*x)-besselj(n-l+1,p0.*x)) ...
                .*(l.*besselj(l,p1.*x)./(p1.*x)-besselj(l+1,p1.*x));
            i16=integral(f16,0,inf);
    %Part 1
            d1=E0x.*n.*(n-l).*wce.^2./nper./n0x./ve.^2;
            d2=-1i.*E0y.*n.*wce./nper./ve;
            d3=-1i.*E0x.*(n-l).*wce./n0x./ve;
            d4=-E0y;
            d7=E0x.*(n-l).*wce./n0x./ve;
            d8=-1i.*E0y;
            b1=-(n0x.*n1y+2i.*(n-l).*wce.^2./ve.^2);
            b2=-(n0x.*n1y+2i.*wce./ve.^2);
            c1=b1.*l.*wce;
            c2=(n-l).*n1per.^2.*wce;
            c3=-2i.*l.*wce.^2./ve;
            c4=n1per.^2.*ve;
            c5=2i.*l.*n0z.*wce.^2./ve;
            axx1=((c1.*cos(de1)+c2.*sin(de1)).*Z1+(c3.*cos(de1)+c4.*sin(de1)).*Z2 ...
                +c5.*cos(de1).*Z3)./n1per.*(d1.*cos(de).*i1+d2.*cos(de).*i2+d3.*sin(de).*i3+d4.*sin(de).*i4);
            axy1=((c1.*sin(de1)-c2.*cos(de1)).*Z1+(c3.*sin(de1)-c4.*cos(de1)).*Z2 ...
                +c5.*sin(de1).*Z3)./n1per.*(d1.*cos(de).*i1+d2.*cos(de).*i2+d3.*sin(de).*i3+d4.*sin(de).*i4);
            axz1=(b1.*Z4-2i.*wce./ve.*Z12+2i.*n0z.*wce./ve.*(-Z13+zl1./2.*Z5)) ...
                .*(d1.*cos(de).*i1+d2.*cos(de).*i2+d3.*sin(de).*i3+d4.*sin(de).*i4);
            ayx1=((c1.*cos(de1)+c2.*sin(de1)).*Z1+(c3.*cos(de1)+c4.*sin(de1)).*Z2 ...
                +c5.*cos(de1).*Z3)./n1per.*(d1.*sin(de).*i1+d2.*sin(de).*i2-d3.*cos(de).*i3-d4.*cos(de).*i4);
            ayy1=((c1.*sin(de1)-c2.*cos(de1)).*Z1+(c3.*sin(de1)-c4.*cos(de1)).*Z2 ...
                +c5.*sin(de1).*Z3)./n1per.*(d1.*sin(de).*i1+d2.*sin(de).*i2-d3.*cos(de).*i3-d4.*cos(de).*i4);
            ayz1=(b1.*Z4-2i.*wce./ve.*Z12+2i.*n0z.*wce./ve.*(-Z13+zl1./2.*Z5)) ...
                .*(d1.*sin(de).*i1+d2.*sin(de).*i2-d3.*cos(de).*i3-d4.*cos(de).*i4);
            azx1=((c1.*cos(de1)+c2.*sin(de1)).*Z4+(c3.*cos(de1)+c4.*sin(de1)).*Z12 ...
                +c5.*cos(de1).*Z14)./n1per.*(d7.*i1+d8.*i2);
            azy1=((c1.*sin(de1)-c2.*cos(de1)).*Z4+(c3.*sin(de1)-c4.*cos(de1)).*Z12 ...
                +c5.*sin(de1).*Z14)./n1per.*(d7.*i1+d8.*i2);
            azz1=(b2.*Z13+1i.*n0z.*wce./ve.*zl1.*Z15).*(d7.*i1+d8.*i2);
            %Next
            c1=1i.*b1.*ve;
            c2=2.*wce;
            c3=-2.*n0z.*wce;
            axx2=sin(de1).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*cos(de).*i5+d2.*cos(de).*i6+d3.*sin(de).*i7+d4.*sin(de).*i8);
            axy2=-cos(de1).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*cos(de).*i5+d2.*cos(de).*i6+d3.*sin(de).*i7+d4.*sin(de).*i8);
            ayx2=sin(de1).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*sin(de).*i5+d2.*sin(de).*i6-d3.*cos(de).*i7-d4.*cos(de).*i8);
            ayy2=-cos(de1).*(c1.*Z1+c2.*Z2+c3.*Z3).*(d1.*sin(de).*i5+d2.*sin(de).*i6-d3.*cos(de).*i7-d4.*cos(de).*i8);
            azx2=sin(de1).*(c1.*Z4+c2.*Z12+c3.*Z14).*(d7.*i5+d8.*i6);
            azy2=-cos(de1).*(c1.*Z4+c2.*Z12+c3.*Z14).*(d7.*i5+d8.*i6);
            %Next
            c1=-l.^2.*wce.^2./n1per./ve;
            c2=(n-l).*wce./ve;
            axx3=c1.*sin(de1).*(Z2+c2.*Z1).*(d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12);
            axy3=-c1.*cos(de1).*(Z2+c2.*Z1).*(d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12);
            ayx3=c1.*sin(de1).*(Z2+c2.*Z1).*(d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12);
            ayy3=-c1.*cos(de1).*(Z2+c2.*Z1).*(d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12);
            azx3=c1.*sin(de1).*(Z12+c2.*Z4).*(d7.*i9+d8.*i10);
            azy3=-c1.*cos(de1).*(Z12+c2.*Z4).*(d7.*i9+d8.*i10);
            %Next
            c1=1i.*l.*wce;
            axx4=c1.*cos(de1).*(Z2+c2.*Z1).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);    %此项为零？
            axy4=c1.*sin(de1).*(Z2+c2.*Z1).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);
            axz4=1i.*n1per.*(Z12+(n-l).*wce./ve.*Z4).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);
            ayx4=c1.*cos(de1).*(Z2+c2.*Z1).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
            ayy4=c1.*sin(de1).*(Z2+c2.*Z1).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
            ayz4=1i.*n1per.*(Z12+(n-l).*wce./ve.*Z4).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
            azx4=c1.*cos(de1).*(Z12+c2.*Z4).*(d7.*i13+d8.*i14);
            azy4=c1.*sin(de1).*(Z12+c2.*Z4).*(d7.*i13+d8.*i14);
            azz4=1i.*n1per.*(Z13./ve-n0z.*Z18).*(d7.*i13+d8.*i14);
            %Sum
            sxx1=30./w0./nz./n1z.*(axx1+axx2+axx3+axx4);
            sxy1=30./w0./nz./n1z.*(axy1+axy2+axy3+axy4);
            sxz1=30./w0./nz./n1z.*ve.*(axz1+axz4);
            syx1=30./w0./nz./n1z.*(ayx1+ayx2+ayx3+ayx4);
            syy1=30./w0./nz./n1z.*(ayy1+ayy2+ayy3+ayy4);
            syz1=30./w0./nz./n1z.*ve.*(ayz1+ayz4);
            szx1=30./w0./nz./n1z.*(azx1+azx2+azx3+azx4);
            szy1=30./w0./nz./n1z.*(azy1+azy2+azy3+azy4);
            szz1=30./w0./nz./n1z.*ve.*(azz1+azz4);
    %Part 2
            d5=n.*wce./nper./ve;
            d6=-1i;
            c1=b2.*l.*wce;
            c2=(n-l).*n1per.^2.*wce;
            c3=1i.*l.*(1-(n-l).*wce).*wce.^2./ve.^2;
            axx1=((c1.*cos(de1)+c2.*sin(de1)).*Z4+c3.*cos(de1).*Z5)./n1per.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            axy1=((c1.*sin(de1)-c2.*cos(de1)).*Z4+c3.*sin(de1).*Z5)./n1per.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            axz1=(b2.*ve.*Z13+1i.*(1-(n-l).*wce).*wce./ve.*zl1.*Z5).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            ayx1=((c1.*cos(de1)+c2.*sin(de1)).*Z4+c3.*cos(de1).*Z5)./n1per.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            ayy1=((c1.*sin(de1)-c2.*cos(de1)).*Z4+c3.*sin(de1).*Z5)./n1per.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            ayz1=(b2.*ve.*Z13+1i.*(1-(n-l).*wce).*wce./ve.*zl1.*Z5).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            azx1=((c1.*cos(de1)+c2.*sin(de1)).*Z13+c3.*cos(de1).*Z15)./n1per.*i1;
            azy1=((c1.*sin(de1)-c2.*cos(de1)).*Z13+c3.*sin(de1).*Z15)./n1per.*i1;
            azz1=(b2.*ve.*Z18+1i.*wce.*(1-(n-l).*wce)./ve.*zl1.*Z15).*i1;
            %Next
            axx2=ve.*sin(de1).*(1i.*b2.*Z4-(1-(n-l).*wce).*wce./ve.^2.*Z5).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
            axy2=-ve.*cos(de1).*(1i.*b2.*Z4-(1-(n-l).*wce).*wce./ve.^2.*Z5).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
            ayx2=ve.*sin(de1).*(1i.*b2.*Z4-(1-(n-l).*wce).*wce./ve.^2.*Z5).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
            ayy2=-ve.*cos(de1).*(1i.*b2.*Z4-(1-(n-l).*wce).*wce./ve.^2.*Z5).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
            azx2=ve.*sin(de1).*(1i.*b2.*Z13-(1-(n-l).*wce).*wce./ve.^2.*Z15).*i5;
            azy2=-ve.*cos(de1).*(1i.*b2.*Z13-(1-(n-l).*wce).*wce./ve.^2.*Z15).*i5;
            %Next
            axx3=-l.^2.*(n-l).*wce.^3./n1per./ve.^2.*sin(de1).*Z4.*(d5.*cos(de).*i9+d6.*sin(de).*i11);
            axy3=l.^2.*(n-l).*wce.^3./n1per./ve.^2.*cos(de1).*Z4.*(d5.*cos(de).*i9+d6.*sin(de).*i11);
            ayx3=-l.^2.*(n-l).*wce.^3./n1per./ve.^2.*sin(de1).*Z4.*(d5.*sin(de).*i9-d6.*cos(de).*i11);
            ayy3=l.^2.*(n-l).*wce.^3./n1per./ve.^2.*cos(de1).*Z4.*(d5.*sin(de).*i9-d6.*cos(de).*i11);
            azx3=-l.^2.*(n-l).*wce.^3./n1per./ve.^2.*sin(de1).*Z13.*i9;
            azy3=l.^2.*(n-l).*wce.^3./n1per./ve.^2.*cos(de1).*Z13.*i9;
            %Next
            axx4=1i.*l.*(n-l).*wce.^2./ve.*cos(de1).*Z4.*(d5.*cos(de).*i13+d6.*sin(de).*i15);    %此项为零?
            axy4=1i.*l.*(n-l).*wce.^2./ve.*sin(de1).*Z4.*(d5.*cos(de).*i13+d6.*sin(de).*i15);
            axz4=1i.*(n-l).*n1per.*wce.*Z13.*(d5.*cos(de).*i13+d6.*sin(de).*i15);
            ayx4=1i.*l.*(n-l).*wce.^2./ve.*cos(de1).*Z4.*(d5.*sin(de).*i13-d6.*cos(de).*i15);
            ayy4=1i.*l.*(n-l).*wce.^2./ve.*sin(de1).*Z4.*(d5.*sin(de).*i13-d6.*cos(de).*i15);
            ayz4=1i.*(n-l).*n1per.*wce.*Z13.*(d5.*sin(de).*i13-d6.*cos(de).*i15);
            azx4=1i.*l.*(n-l).*wce.^2./ve.*cos(de1).*Z13.*i13;
            azy4=1i.*l.*(n-l).*wce.^2./ve.*sin(de1).*Z13.*i13;
            azz4=1i.*(n-l).*n1per.*wce.*Z18.*i13;
            %Sum
            sxx2=30./w0./nz./n1z.*E0z.*(axx1+axx2+axx3+axx4);
            sxy2=30./w0./nz./n1z.*E0z.*(axy1+axy2+axy3+axy4);
            sxz2=30./w0./nz./n1z.*E0z.*(axz1+axz4);
            syx2=30./w0./nz./n1z.*E0z.*(ayx1+ayx2+ayx3+ayx4);
            syy2=30./w0./nz./n1z.*E0z.*(ayy1+ayy2+ayy3+ayy4);
            syz2=30./w0./nz./n1z.*E0z.*(ayz1+ayz4);
            szx2=30./w0./nz./n1z.*E0z.*(azx1+azx2+azx3+azx4);
            szy2=30./w0./nz./n1z.*E0z.*(azy1+azy2+azy3+azy4);
            szz2=30./w0./nz./n1z.*E0z.*(azz1+azz4);
    %Part 3
            d9=l.*wce./n1per./ve;
            axx1=l.*wce./n1per.*cos(de1).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            axy1=l.*wce./n1per.*sin(de1).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            ayx1=l.*wce./n1per.*cos(de1).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            ayy1=l.*wce./n1per.*sin(de1).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            axx2=1i.*ve.*sin(de1).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
            axy2=-1i.*ve.*cos(de1).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
            ayx2=1i.*ve.*sin(de1).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
            ayy2=-1i.*ve.*cos(de1).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
            %Sum
            sxx3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*Z2.*(axx1+axx2);
            sxy3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*Z2.*(axy1+axy2);
            sxz3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*ve.*Z12.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            syx3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*Z2.*(ayx1+ayx2);
            syy3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*Z2.*(ayy1+ayy2);
            syz3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*ve.*Z12.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            szx3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*ve.*Z12.*(d9.*cos(de1).*i1+1i.*sin(de1).*i5);
            szy3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*ve.*Z12.*(d9.*sin(de1).*i1-1i.*cos(de1).*i5);
            szz3=30./w0./nz./n1z.*(n1x.*E0y-n1y.*E0x).*n0z.*ve.*(zn0.*Z13-Z18).*i1;
            %Next
            d1=-1i.*E0y.*n.*(n-l).*wce.^2./nper./n0x./ve.^2;
            d2=E0x.*n.*wce./nper./ve;
            d3=-E0y.*(n-l).*wce./n0x./ve;
            d4=-1i.*E0x;
            d7=-1i.*E0y.*(n-l).*wce./n0x./ve;
            axx1=l.*wce./n1per./ve.*cos(de1).*(d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12);
            axy1=l.*wce./n1per./ve.*sin(de1).*(d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12);
            axz1=d1.*cos(de).*i9+d2.*cos(de).*i10+d3.*sin(de).*i11+d4.*sin(de).*i12;
            ayx1=l.*wce./n1per./ve.*cos(de1).*(d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12);
            ayy1=l.*wce./n1per./ve.*sin(de1).*(d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12);
            ayz1=d1.*sin(de).*i9+d2.*sin(de).*i10-d3.*cos(de).*i11-d4.*cos(de).*i12;
            azx1=d9.*cos(de1).*(d7.*i9+E0x.*i10);
            azy1=d9.*sin(de1).*(d7.*i9+E0x.*i10);
            azz1=d7.*i9+E0x.*i10;
            axx2=1i.*sin(de1).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);
            axy2=-1i.*cos(de1).*(d1.*cos(de).*i13+d2.*cos(de).*i14+d3.*sin(de).*i15+d4.*sin(de).*i16);
            ayx2=1i.*sin(de1).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
            ayy2=-1i.*cos(de1).*(d1.*sin(de).*i13+d2.*sin(de).*i14-d3.*cos(de).*i15-d4.*cos(de).*i16);
            azx2=1i.*sin(de1).*(d7.*i13+E0x.*i14);
            azy2=-1i.*cos(de1).*(d7.*i13+E0x.*i14);
            %Sum
            sxx4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z2+(n-l).*wce./ve.*Z1).*(axx1+axx2);
            sxy4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z2+(n-l).*wce./ve.*Z1).*(axy1+axy2);
            sxz4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z12+(n-l).*wce./ve.*Z4).*axz1;
            syx4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z2+(n-l).*wce./ve.*Z1).*(ayx1+ayx2);
            syy4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z2+(n-l).*wce./ve.*Z1).*(ayy1+ayy2);
            syz4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z12+(n-l).*wce./ve.*Z4).*ayz1;
            szx4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z12+(n-l).*wce./ve.*Z4).*(azx1+azx2);
            szy4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z12+(n-l).*wce./ve.*Z4).*(azy1+azy2);
            szz4=-30./w0./nz./n1z.*1i.*l.*wce.*(Z13./ve-n0z.*Z18).*azz1;
            %Next
            axx1=l.*wce./n1per./ve.*cos(de1).*(d5.*cos(de).*i10+d6.*sin(de).*i12);
            axy1=l.*wce./n1per./ve.*sin(de1).*(d5.*cos(de).*i10+d6.*sin(de).*i12);
            ayx1=l.*wce./n1per./ve.*cos(de1).*(d5.*sin(de).*i10-d6.*cos(de).*i12);
            ayy1=l.*wce./n1per./ve.*sin(de1).*(d5.*sin(de).*i10-d6.*cos(de).*i12);
            axx2=1i.*sin(de1).*(d5.*cos(de).*i14+d6.*sin(de).*i16);
            axy2=-1i.*cos(de1).*(d5.*cos(de).*i14+d6.*sin(de).*i16);
            ayx2=1i.*sin(de1).*(d5.*sin(de).*i14-d6.*cos(de).*i16);
            ayy2=-1i.*cos(de1).*(d5.*sin(de).*i14-d6.*cos(de).*i16);
            %Sum
            sxx5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z4.*(axx1+axx2);
            sxy5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z4.*(axy1+axy2);
            sxz5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z13.*(d5.*cos(de).*i10+d6.*sin(de).*i12);
            syx5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z4.*(ayx1+ayx2);
            syy5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z4.*(ayy1+ayy2);
            syz5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z13.*(d5.*sin(de).*i10-d6.*cos(de).*i12);
            szx5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z13.*(d9.*cos(de1).*i10+1i.*sin(de1).*i14);
            szy5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z13.*(d9.*sin(de1).*i10-1i.*cos(de1).*i14);
            szz5=-30./w0./nz./n1z.*1i.*l.*n0x.*wce.*E0z.*Z18.*i10;
            %Next
            axx1=l.*wce./n1per./ve.*cos(de1).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            axy1=l.*wce./n1per./ve.*sin(de1).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            ayx1=l.*wce./n1per./ve.*cos(de1).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            ayy1=l.*wce./n1per./ve.*sin(de1).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            axx2=1i.*sin(de1).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
            axy2=-1i.*cos(de1).*(d5.*cos(de).*i5+d6.*sin(de).*i7);
            ayx2=1i.*sin(de1).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
            ayy2=-1i.*cos(de1).*(d5.*sin(de).*i5-d6.*cos(de).*i7);
            %Sum
            sxx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z1.*(axx1+axx2);
            sxy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z1.*(axy1+axy2);
            sxz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z4.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            syx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z1.*(ayx1+ayx2);
            syy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z1.*(ayy1+ayy2);
            syz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z4.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            szx6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z4.*(d9.*cos(de1).*i1+1i.*sin(de1).*i5);
            szy6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z4.*(d9.*sin(de1).*i1-1i.*cos(de1).*i5);
            szz6=30./w0./nz./n1z.*l.*n0x.*wce.*E0y.*Z13.*i1;
    %Part 4
            d1=n.*l.*wce.^2./nper./n1per./ve.^2;
            d2=1i.*n.*wce./nper./ve;
            d3=-1i.*l.*wce./n1per./ve;
            d4=1;
            b3=-2i.*(n-l).*wce.^2./n0x./ve.^2.*E0x-n0x.*E0y;
            c1=(n0x.*n1y-2i.*l.*wce.^2./ve.^2).*E0z-2i.*(n-l).*n1z.*wce.^2./n0x./ve.^2.*E0x;
            c2=((n-l).*n1y.*E0x+l.*b3).*wce./ve;
            c3=b3;
            c4=-2i.*wce./ve.*E0z;
            c5=-2i.*n1z.*wce./ve.*E0z;
            c6=1i.*wce.*((1-(n-l).*wce).*E0z+(n-l).*n0z./n0x.*wce.*E0x).*n1z./n0z./ve.^2;
            c7=c1.*Z9+c2.*Z6+c3.*Z7+c4.*Z10+c5.*Z11+c6.*Z8;
            c8=(n0x.*n1y-2i.*wce.*w1./ve.^2).*ve;
            c9=-l.*n0x.*wce.*E0y;
            c10=1i.*wce.*(w1-l.*wce).*((1-(n-l).*wce).*E0z+(n-l).*n0z.*wce./n0x.*E0x)./n0z./ve.^2;
            c11=c8.*(E0z.*Z11+E0x.*(n-l).*wce./n0x./ve.*Z9)+c9.*Z9+c10.*Z8;
            c12=(n0x.*n1y-2i.*w1.*wce./ve.^2).*E0z-2i.*(n-l).*n1z.*wce.^2./n0x./ve.^2.*E0x;
            axx1=c7.*(d1.*cos(de).*cos(de1).*i1+d2.*cos(de).*sin(de1).*i5+d3.*sin(de).*cos(de1).*i3+d4.*sin(de).*sin(de1).*i7);
            axy1=c7.*(d1.*cos(de).*sin(de1).*i1-d2.*cos(de).*cos(de1).*i5+d3.*sin(de).*sin(de1).*i3-d4.*sin(de).*cos(de1).*i7);
            axz1=c11.*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            ayx1=c7.*(d1.*sin(de).*cos(de1).*i1+d2.*sin(de).*sin(de1).*i5-d3.*cos(de).*cos(de1).*i3-d4.*cos(de).*sin(de1).*i7);
            ayy1=c7.*(d1.*sin(de).*sin(de1).*i1-d2.*sin(de).*cos(de1).*i5-d3.*cos(de).*sin(de1).*i3+d4.*cos(de).*cos(de1).*i7);
            ayz1=c11.*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            azx1=(c12.*Z11+c2.*Z9+c3.*Z10+c6.*Z16).*(d9.*cos(de1).*i1+1i.*sin(de1).*i5);
            azy1=(c12.*Z11+c2.*Z9+c3.*Z10+c6.*Z16).*(d9.*sin(de1).*i1-1i.*cos(de1).*i5);
            azz1=(c8.*(E0z.*Z17+E0x.*(n-l).*wce./n0x./ve.*Z11)+c9.*Z11+c10.*Z16).*i1;
            %Next
            c1=(-1i.*n0x.*n1y-2.*l.*wce.^2./ve.^2).*E0y;
            c2=-2.*wce./ve.*E0y;
            c3=n1z.*wce./ve.*E0y;
            c4=c1.*Z6+c2.*Z7+c3.*(-2.*Z9+Z8);
            axx2=c4.*(d1.*cos(de).*cos(de1).*i2+d2.*cos(de).*sin(de1).*i6+d3.*sin(de).*cos(de1).*i4+d4.*sin(de).*sin(de1).*i8);
            axy2=c4.*(d1.*cos(de).*sin(de1).*i2-d2.*cos(de).*cos(de1).*i6+d3.*sin(de).*sin(de1).*i4-d4.*sin(de).*cos(de1).*i8);
            axz2=ve.*E0y.*(-1i.*(n0x.*n1y-2i.*wce.*w1./ve.^2).*Z9+wce.*(w1-l.*wce)./ve.^2.*Z8).*(d5.*cos(de).*i2+d6.*sin(de).*i4);
            ayx2=c4.*(d1.*sin(de).*cos(de1).*i2+d2.*sin(de).*sin(de1).*i6-d3.*cos(de).*cos(de1).*i4-d4.*cos(de).*sin(de1).*i8);
            ayy2=c4.*(d1.*sin(de).*sin(de1).*i2-d2.*sin(de).*cos(de1).*i6-d3.*cos(de).*sin(de1).*i4+d4.*cos(de).*cos(de1).*i8);
            ayz2=ve.*E0y.*(-1i.*(n0x.*n1y-2i.*wce.*w1./ve.^2).*Z9+wce.*(w1-l.*wce)./ve.^2.*Z8).*(d5.*sin(de).*i2-d6.*cos(de).*i4);
            azx2=(c1.*Z9+c2.*Z10+c3.*(-2.*Z11+Z16)).*(d9.*cos(de1).*i2+1i.*sin(de1).*i6);
            azy2=(c1.*Z9+c2.*Z10+c3.*(-2.*Z11+Z16)).*(d9.*sin(de1).*i2-1i.*cos(de1).*i6);
            azz2=ve.*E0y.*(-1i.*(n0x.*n1y-2i.*wce.*w1./ve.^2).*Z11+wce.*(w1-l.*wce)./ve.^2.*Z16).*i2;
            %Next
            c1=(n-l).^2.*wce.^2./n0x./ve.^2.*E0y.*(Z7+l.*wce./ve.*Z6);
            axx3=c1.*(d1.*cos(de).*cos(de1).*i9+d2.*cos(de).*sin(de1).*i13+d3.*sin(de).*cos(de1).*i11+d4.*sin(de).*sin(de1).*i15);
            axy3=c1.*(d1.*cos(de).*sin(de1).*i9-d2.*cos(de).*cos(de1).*i13+d3.*sin(de).*sin(de1).*i11-d4.*sin(de).*cos(de1).*i15);
            axz3=E0y.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z9.*(d5.*cos(de).*i9+d6.*sin(de).*i11);
            ayx3=c1.*(d1.*sin(de).*cos(de1).*i9+d2.*sin(de).*sin(de1).*i13-d3.*cos(de).*cos(de1).*i11-d4.*cos(de).*sin(de1).*i15);
            ayy3=c1.*(d1.*sin(de).*sin(de1).*i9-d2.*sin(de).*cos(de1).*i13-d3.*cos(de).*sin(de1).*i11+d4.*cos(de).*cos(de1).*i15);
            ayz3=E0y.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z9.*(d5.*sin(de).*i9-d6.*cos(de).*i11);
            azx3=(n-l).^2.*wce.^2./n0x./ve.^2.*E0y.*(Z10+l.*wce./ve.*Z9).*(d9.*cos(de1).*i9+1i.*sin(de1).*i13);
            azy3=(n-l).^2.*wce.^2./n0x./ve.^2.*E0y.*(Z10+l.*wce./ve.*Z9).*(d9.*sin(de1).*i9-1i.*cos(de1).*i13);
            azz3=E0y.*l.*(n-l).^2.*wce.^3./n0x./ve.^2.*Z11.*i9;
            %Next
            c1=1i.*n0x.*E0z;
            c2=1i.*l.*n0x.*wce./ve.*E0z;
            c3=1i.*(n-l).*wce./ve.*E0x;
            c4=1i.*l.*(n-l).*wce.^2./ve.^2.*E0x;
            c5=c1.*Z10+c2.*Z9+c3.*Z7+c4.*Z6;
            c6=1i.*n0x.*w1./ve.*E0z;
            c7=-1i.*n0x.*n1z.*E0z;
            axx4=c5.*(d1.*cos(de).*cos(de1).*i10+d2.*cos(de).*sin(de1).*i14+d3.*sin(de).*cos(de1).*i12+d4.*sin(de).*sin(de1).*i16);
            axy4=c5.*(d1.*cos(de).*sin(de1).*i10-d2.*cos(de).*cos(de1).*i14+d3.*sin(de).*sin(de1).*i12-d4.*sin(de).*cos(de1).*i16);
            axz4=1i.*l.*wce.*(n0x.*E0z.*Z11+(n-l).*wce./ve.*E0x.*Z9).*(d5.*cos(de).*i10+d6.*sin(de).*i12);
            ayx4=c5.*(d1.*sin(de).*cos(de1).*i10+d2.*sin(de).*sin(de1).*i14-d3.*cos(de).*cos(de1).*i12-d4.*cos(de).*sin(de1).*i16);
            ayy4=c5.*(d1.*sin(de).*sin(de1).*i10-d2.*sin(de).*cos(de1).*i14-d3.*cos(de).*sin(de1).*i12+d4.*cos(de).*cos(de1).*i16);
            ayz4=1i.*l.*wce.*(n0x.*E0z.*Z11+(n-l).*wce./ve.*E0x.*Z9).*(d5.*sin(de).*i10-d6.*cos(de).*i12);
            azx4=(c6.*Z11+c7.*Z17+c3.*Z10+c4.*Z9).*(d9.*cos(de1).*i10+1i.*sin(de1).*i14);
            azy4=(c6.*Z11+c7.*Z17+c3.*Z10+c4.*Z9).*(d9.*sin(de1).*i10-1i.*cos(de1).*i14);
            azz4=1i.*l.*n0x.*wce.*(E0z.*Z17+E0x.*(n-l).*wce./n0x./ve.*Z11).*i10;
            %Sum
            sxx7=30./w0./nz./n0z.*ve./w1.*(axx1+axx2+axx3+axx4);
            sxy7=30./w0./nz./n0z.*ve./w1.*(axy1+axy2+axy3+axy4);
            sxz7=30./w0./nz./n0z./w1.*(axz1+axz2+axz3+axz4);
            syx7=30./w0./nz./n0z.*ve./w1.*(ayx1+ayx2+ayx3+ayx4);
            syy7=30./w0./nz./n0z.*ve./w1.*(ayy1+ayy2+ayy3+ayy4);
            syz7=30./w0./nz./n0z./w1.*(ayz1+ayz2+ayz3+ayz4);
            szx7=30./w0./nz./n0z.*ve./w1.*(azx1+azx2+azx3+azx4);
            szy7=30./w0./nz./n0z.*ve./w1.*(azy1+azy2+azy3+azy4);
            szz7=30./w0./nz./n0z./w1.*(azz1+azz2+azz3+azz4);
    %Part 5
            axy1=(E0z.*Z10+E0x.*(n-l).*wce./n0x./ve.*Z7).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            ayy1=(E0z.*Z10+E0x.*(n-l).*wce./n0x./ve.*Z7).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            azy1=(E0z.*((w1-l.*wce)./ve.*Z11-n1z.*Z17)+E0x.*(n-l).*wce./n0x./ve.*Z10).*i1;
            axy2=-1i.*E0y.*Z7.*(d5.*cos(de).*i2+d6.*sin(de).*i4);
            ayy2=-1i.*E0y.*Z7.*(d5.*sin(de).*i2-d6.*cos(de).*i4);
            azy2=-1i.*E0y.*Z10.*i2;
            sxy8=30./w0./nz./n0z.*n0x.*ve./w1.*(axy1+axy2);
            syy8=30./w0./nz./n0z.*n0x.*ve./w1.*(ayy1+ayy2);
            szy8=30./w0./nz./n0z.*n0x.*ve./w1.*(azy1+azy2);
            %Next
            c1=E0z.*(Z10+l.*wce./ve.*Z9)+E0x.*(n-l).*wce./n0x./ve.*(Z7+l.*wce./ve.*Z6);
            axx1=c1.*(-d1.*cos(de).*sin(de1).*i9+d2.*cos(de).*cos(de1).*i13-d3.*sin(de).*sin(de1).*i11+d4.*sin(de).*cos(de1).*i15);
            axy1=c1.*(d1.*cos(de).*cos(de1).*i9+d2.*cos(de).*sin(de1).*i13+d3.*sin(de).*cos(de1).*i11+d4.*sin(de).*sin(de1).*i15);
            axz1=(E0z.*Z11+E0x.*(n-l).*wce./n0x./ve.*Z9).*(d5.*cos(de).*i13+d6.*sin(de).*i15);
            ayx1=c1.*(-d1.*sin(de).*sin(de1).*i9+d2.*sin(de).*cos(de1).*i13+d3.*cos(de).*sin(de1).*i11-d4.*cos(de).*cos(de1).*i15);
            ayy1=c1.*(d1.*sin(de).*cos(de1).*i9+d2.*sin(de).*sin(de1).*i13-d3.*cos(de).*cos(de1).*i11-d4.*cos(de).*sin(de1).*i15);
            ayz1=(E0z.*Z11+E0x.*(n-l).*wce./n0x./ve.*Z9).*(d5.*sin(de).*i13-d6.*cos(de).*i15);
            azx1=(E0z.*(w1./ve.*Z11-n1z.*Z17)+E0x.*(n-l).*wce./n0x./ve.*(Z10+l.*wce./ve.*Z9)).*(-d9.*sin(de1).*i9+1i.*cos(de1).*i13);
            azy1=(E0z.*(w1./ve.*Z11-n1z.*Z17)+E0x.*(n-l).*wce./n0x./ve.*(Z10+l.*wce./ve.*Z9)).*(d9.*cos(de1).*i9+1i.*sin(de1).*i13);
            azz1=(E0z.*Z17+E0x.*(n-l).*wce./n0x./ve.*Z11).*i13;
            c1=-1i.*E0y.*(Z7+l.*wce./ve.*Z6);
            axx2=c1.*(-d1.*cos(de).*sin(de1).*i10+d2.*cos(de).*cos(de1).*i14-d3.*sin(de).*sin(de1).*i12+d4.*sin(de).*cos(de1).*i16);
            axy2=c1.*(d1.*cos(de).*cos(de1).*i10+d2.*cos(de).*sin(de1).*i14+d3.*sin(de).*cos(de1).*i12+d4.*sin(de).*sin(de1).*i16);
            axz2=-1i.*E0y.*Z9.*(d5.*cos(de).*i14+d6.*sin(de).*i16);
            ayx2=c1.*(-d1.*sin(de).*sin(de1).*i10+d2.*sin(de).*cos(de1).*i14+d3.*cos(de).*sin(de1).*i12-d4.*cos(de).*cos(de1).*i16);
            ayy2=c1.*(d1.*sin(de).*cos(de1).*i10+d2.*sin(de).*sin(de1).*i14-d3.*cos(de).*cos(de1).*i12-d4.*cos(de).*sin(de1).*i16);
            ayz2=-1i.*E0y.*Z9.*(d5.*sin(de).*i14-d6.*cos(de).*i16);
            azx2=-1i.*E0y.*(Z10+l.*wce./ve.*Z9).*(-d9.*sin(de1).*i10+1i.*cos(de1).*i14);
            azy2=-1i.*E0y.*(Z10+l.*wce./ve.*Z9).*(d9.*cos(de1).*i10+1i.*sin(de1).*i14);
            azz2=-1i.*E0y.*Z11.*i14;
            %Sum
            sxx9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(axx1+axx2);
            sxy9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(axy1+axy2);
            sxz9=-30./w0./nz./n0z.*1i.*(n-l).*n1per.*wce./w1.*(axz1+axz2);
            syx9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(ayx1+ayx2);
            syy9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(ayy1+ayy2);
            syz9=-30./w0./nz./n0z.*1i.*(n-l).*n1per.*wce./w1.*(ayz1+ayz2);
            szx9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(azx1+azx2);
            szy9=-30./w0./nz./n0z.*(n-l).*wce./w1.*(azy1+azy2);
            szz9=-30./w0./nz./n0z.*1i.*(n-l).*n1per.*wce./w1.*(azz1+azz2);
            %Next
            a1=(E0z.*Z9+E0x.*(n-l).*wce./n0x./ve.*Z6).*(d5.*cos(de).*i1+d6.*sin(de).*i3);
            ay1=(E0z.*Z9+E0x.*(n-l).*wce./n0x./ve.*Z6).*(d5.*sin(de).*i1-d6.*cos(de).*i3);
            az1=(E0z.*Z11+E0x.*(n-l).*wce./n0x./ve.*Z9).*i1;
            a2=-1i.*E0y.*Z6.*(d5.*cos(de).*i2+d6.*sin(de).*i4);
            ay2=-1i.*E0y.*Z6.*(d5.*sin(de).*i2-d6.*cos(de).*i4);
            az2=-1i.*E0y.*Z9.*i2;
            %Sum
            sxx10=-30./w0./nz./n0z.*(n-l).*n1y.*wce./w1.*(a1+a2);
            sxy10=30./w0./nz./n0z.*(n-l).*n1x.*wce./w1.*(a1+a2);
            syx10=-30./w0./nz./n0z.*(n-l).*n1y.*wce./w1.*(ay1+ay2);
            syy10=30./w0./nz./n0z.*(n-l).*n1x.*wce./w1.*(ay1+ay2);
            szx10=-30./w0./nz./n0z.*(n-l).*n1y.*wce./w1.*(az1+az2);
            szy10=30./w0./nz./n0z.*(n-l).*n1x.*wce./w1.*(az1+az2);
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
            SigNLxx=SigNLxx+exp(1i.*(n.*de-l.*de1)).*sxx;
            SigNLxy=SigNLxy+exp(1i.*(n.*de-l.*de1)).*sxy;
            SigNLxz=SigNLxz+exp(1i.*(n.*de-l.*de1)).*sxz;
            SigNLyx=SigNLyx+exp(1i.*(n.*de-l.*de1)).*syx;
            SigNLyy=SigNLyy+exp(1i.*(n.*de-l.*de1)).*syy;
            SigNLyz=SigNLyz+exp(1i.*(n.*de-l.*de1)).*syz;
            SigNLzx=SigNLzx+exp(1i.*(n.*de-l.*de1)).*szx;
            SigNLzy=SigNLzy+exp(1i.*(n.*de-l.*de1)).*szy;
            SigNLzz=SigNLzz+exp(1i.*(n.*de-l.*de1)).*szz;
        end
    end
    %{
    y=[SigNLxx SigNLxy SigNLxz
        SigNLyx SigNLyy SigNLyz
        SigNLzx SigNLzy SigNLzz];
    %}
    y = 0.096./81.9.*wpe.^2./w./wce./ve.^2.*[SigNLxx SigNLxy SigNLxz
        SigNLyx SigNLyy SigNLyz
        SigNLzx SigNLzy SigNLzz];
    num = num+1;
end
y = y;