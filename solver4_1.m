function [X,Y,Q] = solver4_1(m,n,nt,s)
    %Solver of heat equation in 2D for structured mesh with variable
    % coefficients.
    %   Input:
    %       m: number of elements in x direction
    %       n: number of elements in y direction
    %       nt: number of time steps
    %       s: Source term (1 or 2)
    %   Output:
    %       X,Y: Mesh grid
    %       Q: Solution in matrix form (m*n,nt)
    Lx = 1;
    Ly = 1;
    dx = Lx/m;
    dy = Ly/n;
    tf = 1;
    dt = tf/nt;

    xm = linspace(dx/2,Lx-dx/2,m);
    ym = linspace(dy/2,Ly-dy/2,n);

    [X,Y] = meshgrid(xm,ym);
    if s==1
        S1 = @(x,y) exp(-((x-0.5)^2+(y-0.5)^2)/0.2^2);
        S = arrayfun(S1,X,Y);
        S = reshape (S,m*n,1);
   
    else 
        S2 = source2(X,Y);
        S = 2*reshape(S2,m*n,1);
        S3 = zeros(m*n,1);
       
    end

    In = sparse(eye(n));
    Im = sparse(eye(m));
    
    %example for smooth and positive function a(y)
    a = @(y) y+0.5;
    af = arrayfun(a,Y(:,1));
    A = sparse(diag(af));
    %example for smooth and positive function b(x)
    b = @(x) x^2+.4;
    bf = arrayfun(b,X(1,:));
    B = sparse(diag(bf));
    
    e = ones(m,1);
    Tx = spdiags([e -2*e e],[-1 0 1],m,m);
    Tx(1,1) = -1;
    Tx(m,m) = -1;
    Tx = Tx/dx^2;

    e = ones(n,1);
    Ty = spdiags([e -2*e e],[-1 0 1],n,n);
    Ty(1,1) = -1;
    Ty(n,n) = -1;
    Ty = Ty/dy^2;

    M = -dt*(kron(A',Tx) + kron(Ty',B)) + kron(In',Im);
    [L,U,P] = lu(M);

    Q = zeros(n*m,nt);
    b = zeros(n*m,1);

    for t = 2:nt
        if (s==2) && (t*dt>=0.25)
            S = S3;
        end
        b = Q(:,t-1) + dt*S;
        Q(:,t) = U\(L\(P*b));
    end
end

