%%%Z(zeta)近似解,迭代精度一定，epsilon = 1e-9;包括连分式的近似

function out=Z(zeta)
epsilon = 1e-9;
% epsilon = 1e-15;
n = 0;
if(abs(imag(zeta))<1)
    y = 1i.*sqrt(pi).*exp(-zeta.^2);
    yold = y-10;
%     fprintf('ystart= %d',y); 

%     if(abs(real(zeta))<4)
    if(abs(real(zeta))<4.5)
        while abs(y-yold)>epsilon
       yold = y;
       y =  y + (- 2.* zeta.*(-2.*zeta.^2).^n./(prod(1:2:(2*n+1))));
       n = n+1;
        end 
   
%     elseif(abs(real(zeta))>=4)  %此时虚部肯定小于实部
     elseif(abs(real(zeta))>=4.5)  %此时虚部肯定小于实部       
%         fprintf('imag(zeta)= %f',imag(zeta));  
        while abs(y-yold)>epsilon  %
        yold = y;
%         y = y + (- 1./zeta.*abs(2.*n-1)./(2.*zeta.^2).^n);
        y = y + (- 1./zeta.*(prod(1:2:(abs(2.*n-1))))./(2.*zeta.^2).^n);   %修正
        n = n+1;     
        end 

    end

else 
    zetanew = real(zeta)+1i.*abs(imag(zeta));
    a = 1:11;
    b = 1:11;
    A = 1:11;
    B = 1:11;
    a(1)=zetanew;
    b(1)=-zetanew.^2+0.5;
    for ii=2:11
        a(ii)=-(ii-1).*(2.*(ii-1)-1)./2;     %正负修正！中文的文献中有误！
        b(ii)=-zetanew.^2+0.5+2.*(ii-1);
%         fprintf('bim(ii)= %f',imag(b(ii)));
    end
    A0 = 0;
    Af1 = 1;
    Bf1 = 0;
    B0 = 1;

    A(1) = b(1).*A0 + a(1).*Af1;
    B(1) = b(1).*B0 + a(1).*Bf1;
    A(2) = b(2).*A(1) + a(2).*A0;
    B(2) = b(2).*B(1) + a(2).*B0;    
    for ii=3:11
        A(ii)=b(ii).*A(ii-1) + a(ii).*A(ii-2);
        B(ii)=b(ii).*B(ii-1) + a(ii).*B(ii-2);
%           fprintf('Bim(ii)= %f',imag(B(ii)));    %这里会出现NaN的值
%           fprintf('Aim(ii)= %f',imag(A(ii)));    
    end
%     xx=1:110;
%     plot(xx,A,xx,B);
    y = A(11)/B(11);
    if(imag(zeta)<=-1)
        y =conj(y)+2i.*sqrt(pi).*exp(-zeta^2);        %虚部为负的处理
    end
    
end
% if((abs(real(zeta))>4.01)&&(abs(real(zeta))<=4.5))
%        fprintf('zeta= %d',zeta);  
%        if(abs(imag(zeta))>=1)
%        fprintf('zetaIm= %d',imag(zeta));   
%        end
%        fprintf('n= %d',n); 
% end
out = y;

% toc;

