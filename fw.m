function y = fw(wx)
% w = w.*1e-3;
w = (wx(1)+wx(2)*1i).*1e-3;
y = -cond(PIM(w));%the opposite number of the condition number of the dispersion matrix
% % 最小非负特征值
% PIMeig=eig(PIM(w));
% leig=length(PIMeig);
% for ll=1:leig
%     PIMeigM(ll)=abs(PIMeig(ll));
% end
% MinEig=min(PIMeigM);
% for ll=1:leig
%     if(PIMeigM(ll)==MinEig)
%         mm=ll;
%     end
end
% %y=[real(PIMeig(mm));imag(PIMeig(mm))];
% y = PIMeig(mm);

% 行列式
%y = det(PIM(w));

% 标量
%y = PIM(w);