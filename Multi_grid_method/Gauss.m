function u = Gauss(A,u0,f)
% Perform one Gauss-Seidel iteration on Au = f using initial guss u0

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

A_new = D-L;
f_new = U*u0+f;
u = ForwardSubstitution(A_new, f_new);
end
