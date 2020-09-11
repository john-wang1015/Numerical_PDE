function [u,e] = GSSequence(A,u0,f,tol,maxit)
% Perform a sequence of Gauss-Seidel iteration on Au = f
% using the initial estimate u0, until the 2-norm of the residual
% has been reduced by a factor of at least tol, or until the number 
% of iterations reaches maxit

r = f - A*u0;
vold = u0;
i = 0;
err = [];

while (i < maxit) && (norm(r) > tol)
    vnew = Gauss(A,vold,f);
    err = [err norm(vnew-A\f)];
    vold = vnew;
    r = f-A*vold;
    i = i + 1;
end
e = err;
u = vold;

end
