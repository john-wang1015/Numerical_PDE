function err = VerifyWJ(N, omega)

% Generate a random matrix with strong diagonal dominance
A   = rand(N) + N*eye(N);
f   = rand(N,1);

% Perform GS relaxation
% zero initial guess
u0  = zeros(N,1);
% residual reduction tolerance
tol = 1.e-8;
maxit   = 100;

uWJ = WJSequence(A, u0, f, omega, tol, maxit);

% Compare with the built-in MATLAB solver
err = norm(uWJ - A\f);
