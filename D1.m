%下边频波的包含介电张量的线性项
function y = D1(w)
w1 = w-1;
global n1squ n1x n1y n1z;
n1n1Mat = [n1x.^2 n1x.*n1y n1x.*n1z
    n1y.*n1x n1y.^2 n1y.*n1z
    n1z.*n1x n1z.*n1y n1z.^2];
% 电磁
y = feps1(w)+(n1n1Mat-n1squ.*eye(3))./w1.^2;

% 静电
%y = feps1(w);