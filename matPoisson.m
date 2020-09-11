%-------------------------
function A = matPoisson(N)
% Builds the problem matrix A for NxN unknowns
% (interior points) on a grid with grid spacing h=1/(N+1).
h = 1./(N+1);
e = ones(N^2,1)/h^2;
A = spalloc(N^2, N^2, 5*N^2);
A = spdiags([-e -e], [-N N], N^2, N^2);
b = spdiags([-e 4.*e -e],[-1 0 1],N,N);

for i=1:N
    A(1+(i-1)*N:i*N,1+(i-1)*N:i*N)=b;
end

%  full(a)
sparse(A);

%  figure(5)
%  spy(a)
