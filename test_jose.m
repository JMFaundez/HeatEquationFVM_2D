clear all; close all
m = 50;
Nt = 80;
source = 2;

[X,Y,Q] = solver(m,m,Nt,source);

Qs = sum(Q,1);

t = linspace(0,1,Nt);
if source ==1
  Qa = t.*pi*0.2^2*(erf(0.5/0.2))^2;
elseif source ==2
  Qa = 2*t;
  Qa(t>0.25)=2*0.25;
end

Qs = Qs./(m^2);

figure()
semilogy(t,abs(Qa-Qs)./abs(Qa))
xlabel('time')
ylabel('Absolute Error')
hold off

Qm = Q(:,10);
figure()
surf(X,Y,reshape(Qm,m,m))
title('N=10')

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
