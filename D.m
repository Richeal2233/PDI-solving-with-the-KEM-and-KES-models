%��Ƶģ�İ������������������
function y = D(w)
global nsqu nx ny nz;
nnMat = [nx.^2 nx.*ny nx.*nz
    ny.*nx ny.^2 ny.*nz
    nz.*nx nz.*ny nz.^2];
% ���
y = feps(w)+(nnMat-nsqu.*eye(3))./w.^2;

%����
%y = feps(w);