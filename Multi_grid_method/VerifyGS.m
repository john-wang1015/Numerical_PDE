function err = VerifyGS(N)

% Generate a random matrix with strong diagonal dominance
A   = rand(N) + N*eye(N);
f   = rand(N,1);

% Perform GS relaxation
% zero initial guess
u0  = zeros(N,1);
% residual reduction tolerance
tol = 1.e-8;
maxit   = 100;

uGS = GSSequence(A, u0, f, tol, maxit);

% Compare with the built-in MATLAB solver
err = norm(uGS - A\f);