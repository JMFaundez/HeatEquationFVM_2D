function del2 = source2(X,Y)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
dx = abs(X(1,1)-X(1,2));
dy = abs(Y(1,1)-Y(2,1));
eps = sqrt(max(dx,dy));
R = sqrt((X-0.5).^2+(Y-0.5).^2);
del2 = R*0;
del = pi/(eps^2*(pi^2-4))*(1+cos(pi/eps*R));
del2(R<eps) = del(R<eps);
end

