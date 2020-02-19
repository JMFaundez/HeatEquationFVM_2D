clear all
m = 20;
n = 20;
Lx = 1;
Ly = 1;
dx = Lx/m;
dy = Ly/n;
tf = 2;
nt = 100;
dt = tf/nt;

x = linspace(0,Lx,m+1);
y = linspace(0,Ly,n+1);

xm = linspace(dx/2,Lx-dx/2,m);
ym = linspace(dy/2,Ly-dy/2,n);

[X,Y] = meshgrid(xm,ym);

S1 = @(x,y) exp(-((x-0.5)^2+(y-0.5)^2)/0.2^2);
S1a = arrayfun(S1,X,Y);
S1a = reshape (S1a,m*n,1);

In = sparse(eye(n));
Im = sparse(eye(m));
e = ones(m,1);
Tx = spdiags([e 2*e e],[-1 0 1],m,m);
Tx(1,1) = -1;
Tx(m,m) = -1;
Txf = full(Tx);

e = ones(n,1);
Ty = spdiags([e 2*e e],[-1 0 1],n,n);
Ty(1,1) = -1;
Ty(n,n) = -1;
Tyf = full(Ty);

A = -dt*(kron(In',Tx) + kron(Ty',Im)) + kron(In',Im);
[L,U,P] = lu(A);

Q = zeros(n*m,nt);
b = zeros(n*m,1);

for t = 2:nt
    b = Q(:,t-1) + dt*S1a;
    Q(:,t) = U\(L\(P*b));
end

Qr = reshape(Q(:,nt),m,n);
surf(X,Y,Qr);






