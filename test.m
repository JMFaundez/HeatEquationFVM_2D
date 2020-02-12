m = 50;
n = 50;
nt = 100;
fname = '1.mat';
[X,Y,Q] = solver(m,n,nt,2);
save(fname,'X','Y','Q');
File = load(fname);
Qload = File.Q;
Xload = File.X;
Yload = File.Y;
Q1 = reshape(Qload(:,2),m,n);
surf(Xload,Yload,Q1)