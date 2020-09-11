function vnew = wJacobi(A,vold,f,omega)
% Performs one weighted Jacobi iteration on Au = f with weight omega
% vold is the approximation before the Jacobi iteration, and vnew after
% use sparse matrices for all steps and return in matrix form 

N = length(A);
x1 = diag(A);
D = diag(x1);
A = A - D;
U = zeros(N,N);

for i = 1:N-1
    for j = i+1:N
        U(i,j) = A(i,j);
    end
end

L = -1*(A - U);
U = -1*U;
D_inv = diag(1./x1);

Rj = D_inv*(L+U);
Rw = (1-omega)*eye(N)+omega*Rj;

vnew = Rw*vold + omega*D_inv*f;

end
