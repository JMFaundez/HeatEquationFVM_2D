clear all; close all
m = 40;
Nt = 100;
source = 1;

[X,Y,Q] = solver(m,m,Nt,source);
[X,Y,Q2] = solverLU(m,m,Nt,source);

Qs = sum(Q,1);
Qs2 = sum(Q2,1);

t = linspace(0,1,Nt);
if source ==1
  Qa = t.*pi*0.2^2*(erf(0.5/0.2))^2;
elseif source ==2
  Qa = 2*t;
  Qa(t>0.25)=2*0.25;
end

Qs = Qs./(m^2);
Qs2 = Qs2./(m^2);

figure()
semilogy(t,abs(Qa-Qs))
hold on
semilogy(t,abs(Qa-Qs2))
xlabel('time')
ylabel('Absolute Error')
legend('Previous solver','Fixed by Sophia' )
hold off

Qm = Q(:,70);
figure()
surf(X,Y,reshape(Qm,m,m))
title('N=70')

Qm = Q(:,80);
figure()
surf(X,Y,reshape(Qm,m,m))
title('N=80')

Qm = Q(:,90);
figure()
surf(X,Y,reshape(Qm,m,m))
title('N=90')
%Qm2 = Q2(:,Nt/2);
%figure()
%surf(X,Y,reshape(Qm2,m,m))
%title('Fixed by Sophia')
