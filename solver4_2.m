function [X,Y,Q] = solver4_2(m,n,nt,s)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    Lx = 1;
    Ly = 1;
    dx = Lx/m;
    dy = Ly/n;
    tf = 0.5;
    dt = tf/nt;
    xm = linspace(dx/2,Lx-dx/2,m);
    ym = linspace(dy,Ly-dy,n);

    [X,Y] = meshgrid(xm,ym);
    if s==1
        S1 = @(x,y) 10*exp(-((x-0.5)^2+(y-0.5)^2)/0.2^2);
        S = arrayfun(S1,X,Y);
        S = reshape (S',m*n,1);
    else
        S2 = source2(X,Y);
        S = 2*reshape(S2,m*n,1);
        S3 = zeros(m*n,1);
    end

    In = sparse(eye(n));
    Im = sparse(eye(m));
    e = ones(m,1);
    Tx = spdiags([e -2*e e],[-1 0 1],m,m);
    Tx(1,1) = -1;
    Tx(m,m) = -1;
    Tx = Tx/dx^2;
    %Txf = full(Tx);

    C = zeros(m,n);
    C(1,:) = 1;
    C(end,:) = -1;
    C = reshape(C./dx,m*n,1);

    e = ones(n,1);
    Ty = spdiags([e -2*e e],[-1 0 1],n,n);
    Ty(1,1) = -2;
    Ty(n,n) = -2;
    Ty = Ty/dy^2;
    %Tyf = full(Ty)

    D = zeros(m,n);
    q0 =  sin(pi*xm)/pi;
    q1 = sin(3*pi*xm)/(3*pi) + 1;
    D(:,1) = reshape(q0/dy^2,m,1);
    D(:,n) = reshape(q1/dy^2,m,1);
    D = reshape(D,m*n,1);

    G = C + D + S;

    A = -dt*(kron(In',Tx) + kron(Ty',Im)) + kron(In',Im);
    [L,U,P] = lu(A);
    Q = zeros(n*m,nt);
    b = zeros(n*m,1);

    for t = 2:nt
        if (s==2) && (t*dt>=0.25)
          S = S3;
          G = C + D + S;
        end
            
        b = Q(:,t-1) + dt*G;
        Q(:,t) = U\(L\(P*b));
       
    end
end
