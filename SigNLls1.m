function y = SigNLls1(w)
epsilon = 1e-09;
nl = 0;
y = SigNLls1p5(w,nl);
yold = y-100;
nl = 1;
while max(max(abs(y-yold))) >= epsilon
    yold = y;
    y = SigNLls1p5(w,nl);
    nl = nl+1;
end