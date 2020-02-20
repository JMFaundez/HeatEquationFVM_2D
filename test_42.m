
clear all; close all
m = 40;
Nt = 80;
source = 1;

[X,Y,Q] = solver(m,m,Nt,source);

Qs1 = sum(Q,1);

Qm = Q(:,Nt);
figure()
surf(X,Y,reshape(Qm,m,m))
title('Previous Solver')
Qp = Qm;

[X,Y,Q] = solver4_2(m,m,Nt,source);

Qs2 = sum(Q,1);


%t = linspace(0,1,Nt);
%figure()
%hold on
%plot(t,Qs1-Qs2)
%hold off

Qm = Q(:,Nt/10);
%figure()
%surf(X,Y,reshape(100*abs(Qm-Qp)./Qp,m,m))
%xlabel('x')
%ylabel('y' )
%title('Difference')

figure()
surf(X,Y,reshape(Qm,m,m))
title('Solver 4.2')
xlabel('x')
ylabel('y')
Qp = Qm;
