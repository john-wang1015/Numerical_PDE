function [u,e] = WJSequence(A,u0,f,omega,tol,maxit)
% Perform a sequence of weighted Jacobi itrations on 
% Au = f using the initial estimate u0, until the 2-norm
% of the residual has been reduced by a factor of at least tol, or until
% the number of iteration reaches maxit.

r = f - A*u0;
vold = u0;
i = 0;
err = [];

while (i < maxit) && (norm(r) > tol)
    vnew = wJacobi(A,vold,f,omega);
    err = [err norm(vnew-A\f)];
    vold = vnew;
    r = f-A*vold;
    i = i + 1;
end

e = err;
u = vold;

end