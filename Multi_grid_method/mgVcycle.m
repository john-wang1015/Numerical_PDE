function vnew = mgVcycle(vold,f,alpha1,alpha2,omega,Nmin,N,smoother)
% This is V-cycle implementation, the central part of the code.
% You should base your implementation on the recurisive implementation.
%
% In this routine, first do relaxation on the current level,
% which has N*N interior points, and then recurisiverly call 
% the next level, with (N-1)/2 * (N-1)/2 interior points, followed
% by the coarse grid correction and another relaxation.
%
% > vold is the approxiamtion on the current level before the V-cycle,
%   and vnew is the new approxiamtion after the V-cycle.
%
% > f is the right hand side on the current level, on all levels but 
%   the finest, this is the residual vector.
%
% > alpha1 and alpha2 are the number of pre- and post-relaxations.
%
% > omega is the weight for the weighted Jacobi relaxation.
%
% > if N <= Nmin can use v=A\f solve directly.
%
% > smoother = 1 means weighted Jacobi and smoother = 2 means
% Gauss-Seidel.

A = matPoisson(N);

if N <= Nmin
    vnew = A\f;
else
    if smoother == 1
        for i =1:alpha1
            vnew = wJacobi(A,vold,f,omega);
            vold = vnew;
        end

        M = 0.25*interpolate(N,(N-1)/2)';
        r = f-A*vold;
        error = M*r;
        v2 = zeros(size(error,1),1);
        v2kh = mgVcycle(v2,error,alpha1,alpha2,omega,Nmin,(N-1)/2,smoother);
        v_temp = vold + interpolate(N,(N-1)/2)*v2kh;
    
        for i = 1:alpha2
            vnew = wJacobi(A,v_temp,f,omega);
            v_temp = vnew;
        end
        vnew = v_temp;
    elseif smoother == 2
        for i =1:alpha1
            vnew = Gauss(A,vold,f);
            vold = vnew;
        end

        M = 0.25*interpolate(N,(N-1)/2)';
        error = M*(f-A*vold);
        v2 = zeros(size(error,1),1);
        v2kh = mgVcycle(v2,error,alpha1,alpha2,omega,Nmin,(N-1)/2,smoother);
        v_temp = vold + interpolate(N,(N-1)/2)*v2kh;
    
        for i = 1:alpha2
            vnew = Gauss(A,vold,f);
            v_temp = vnew;
        end
        vnew = v_temp;        
    end
end

end
