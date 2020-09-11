function f = fPossion(N)
% Builds the right hand side column vector f for N*N unknowns
% (iterior points) on a grid with grid spacing h = 1/(N + 1)
h = 1/(N+1);
x = 0:h:1;
y = 0:h:1;
temp = [];

for i = 1:N
    t = zeros(N,1);
    for j = 1:N
        t(j,1) = 2*((1-6*x(i)^2)*y(j)^2*(1-y(j)^2) + (1-6*y(j)^2)*x(i)^2*(1-x(i)^2));
    end
    temp = [temp;t];
end

f = temp;

end
