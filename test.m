clear all; close all
m = [21, 51, 101];
Nt = [50,100,200];
Qz = zeros(length(Nt),length(m));

count = 1;
for t=1:length(Nt)
    for i=1:length(m)
    fname = [num2str(count),'S1.mat'];
    File = load(fname);
    Qload = File.Q;
    Xload = File.X;
    Yload = File.Y;
    Qt = sum(Qload,1);
    Q1 = reshape(Qload(:,end),m(i),m(i));
    Qz(t,i) = Q1((m(i)+1)/2, (m(i)+1)/2);
    count = count+1;
    end
end

t = linspace(0,1,200);
Qa = t.*pi*0.2^2*(erf(0.5/0.2))^2;
Qt = Qt./(m(end)^2);
figure()
hold on
plot(t,Qa-Qt)

dr = 1./m;
dt = 1./Nt;

mr = (Qz(:, 1:end-2)-Qz(:, 2:end-1))./(Qz(:,2:end-1)-Qz(:,3:end));
r = log(mr)
mp = (Qz(1:end-2,:)-Qz(2:end-1,:))./(Qz(2:end-1,:)-Qz(3:end,:));
p = log(mp)

