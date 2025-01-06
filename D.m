%低频模的包含介电张量的线性项
function y = D(w)
global nsqu nx ny nz;
nnMat = [nx.^2 nx.*ny nx.*nz
    ny.*nx ny.^2 ny.*nz
    nz.*nx nz.*ny nz.^2];
% 电磁
y = feps(w)+(nnMat-nsqu.*eye(3))./w.^2;

%静电
%y = feps(w);