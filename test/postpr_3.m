clear all; close all
m = [21, 41, 81, 161, 321];
Nt = [50,100,200,400];
Qz = zeros(length(Nt),length(m));
source = 2;
count = 1;
for t=1:length(Nt)
    for i=1:length(m)
    fname = [num2str(count),'S2.mat'];
    File = load(fname);
    Qload = File.Q;
    Xload = File.X;
    Yload = File.Y;
    Qt = sum(Qload,1);
    index = round(0.1*Nt(t));
    Q1 = reshape(Qload(:,index),m(i),m(i));
    Qz(t,i) = Q1(1, 1);
    count = count+1;
    end
end
Qm = Qload(:,300);
figure()
surf(Xload,Yload,reshape(Qm,m(end),m(end)))
title('N=70')

t = linspace(0,1,Nt(end));
if source ==1
  Qa = t.*pi*0.2^2*(erf(0.5/0.2))^2;
elseif source ==2
  Qa = 2*t;
  Qa(t>=0.25)=2*0.25;
end
Qt = Qt./(m(end)^2);
figure()
plot(t,100*abs((Qa-Qt)./Qa),'LineWidth',2)
ylabel('Relative Error %')
xlabel('time')
saveas(gcf,'p3_s2.png')

dr = 1./m;
dt = 1./Nt;

mr = (Qz(:, 1:end-2)-Qz(:, 2:end-1))./(Qz(:,2:end-1)-Qz(:,3:end));
r = log2(abs(mr))
mp = (Qz(1:end-2,:)-Qz(2:end-1,:))./(Qz(2:end-1,:)-Qz(3:end,:));
p = log2(abs(mp))

error = abs(Qz - Qz(end,end));
%figure()
%loglog(1./m, error(1,:))
%hold on
%loglog(1./m, error(2,:))
%loglog(1./m, error(3,:))
