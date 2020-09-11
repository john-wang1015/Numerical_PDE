%function mgMain
% This is the main driver routine of your code.
% You will first set parameters (including omega, alpha1, alpha2, Nmin,
% number of V-cycles iter, smoother, fine level problem size N, etc.)
% You will create the grid points, the exact PDE solution interpolated
% on the grid, and set up the fine level matrix A and right hand side f.
% Build a for loop for the V-cycles, calling mgVcycle iter times.
% At the end, plot, on the unit square grid, the exact solution of the
% PDE, your multigrid solution for Au=f, and the difference between them. 
% (The ?reshape? command may be helpful to go from lexicographic vectors 
% to 2D arrays.)
% During the V-cycle iterations, keep track of the norms of residuals
% and errors, and output a table similar to Table 4.1 of the textbook.
clc;clear;
N = 31;
A = matPoisson(N);
f = fPossion(N);
Nmin = 3;
omega = 2/3;     
alpha1 = 2;
alpha2 = 1;
smoother = 1;
vold = zeros(N^2,1);

h = 1/(N+1);
x = 0:h:1;
y = 0:h:1;
temp = [];
for i = 1:N
    t = zeros(N,1);
    for j = 1:N
        t(j,1) = (x(i)^2-x(i)^4) * (y(j)^4-y(j)^2);
    end
    temp = [temp;t];
end

vexact = temp;

rh = zeros(15,1);
eh = zeros(15,1);

for i = 1:15
   vnew = mgVcycle(vold,f,alpha1,alpha2,omega,Nmin,N,smoother);
   e = vold - vnew;
   vold = vnew;
   r = A*e;
   eh(i) = norm(e);
   rh(i) = norm(r);
end

r_ratio = zeros(15,1);
e_ratio = zeros(15,1);

for i = 2:15
    r_ratio(i) = rh(i)/rh(i-1);
    e_ratio(i) = eh(i)/eh(i-1);
end

table1 = [rh,r_ratio,eh,e_ratio];

T1 = array2table(table1)
%%
clc;clear;
N = 63;
A = matPoisson(N);
f = fPossion(N);
Nmin = 3;
omega = 2/3;     
alpha1 = 2;
alpha2 = 1;
smoother = 1;
vold = zeros(N^2,1);

h = 1/(N+1);
x = 0:h:1;
y = 0:h:1;
temp = [];
for i = 1:N
    t = zeros(N,1);
    for j = 1:N
        t(j,1) = (x(i)^2-x(i)^4) * (y(j)^4-y(j)^2);
    end
    temp = [temp;t];
end

vexact = temp;

rh = zeros(15,1);
eh = zeros(15,1);

for i = 1:15
   vnew = mgVcycle(vold,f,alpha1,alpha2,omega,Nmin,N,smoother);
   e = vold - vnew;
   vold = vnew;
   r = A*e;
   eh(i) = norm(e);
   rh(i) = norm(r);
end

r_ratio = zeros(15,1);
e_ratio = zeros(15,1);

for i = 2:15
    r_ratio(i) = rh(i)/rh(i-1);
    e_ratio(i) = eh(i)/eh(i-1);
end

table2 = [rh,r_ratio,eh,e_ratio];

T2 = array2table(table2)
%%
N = 31;
A = matPoisson(N);
f = fPossion(N);
Nmin = 3;
omega = 2/3;     
alpha1 = 2;
alpha2 = 1;
smoother = 2;
vold = zeros(N^2,1);

h = 1/(N+1);
x = 0:h:1;
y = 0:h:1;
temp = [];
for i = 1:N
    t = zeros(N,1);
    for j = 1:N
        t(j,1) = (x(i)^2-x(i)^4) * (y(j)^4-y(j)^2);
    end
    temp = [temp;t];
end

vexact = temp;

rh = zeros(15,1);
eh = zeros(15,1);

for i = 1:15
   vnew = mgVcycle(vold,f,alpha1,alpha2,omega,Nmin,N,smoother);
   e = vold - vnew;
   vold = vnew;
   r = A*e;
   eh(i) = norm(e);
   rh(i) = norm(r);
end

r_ratio = zeros(15,1);
e_ratio = zeros(15,1);

for i = 2:15
    r_ratio(i) = rh(i)/rh(i-1);
    e_ratio(i) = eh(i)/eh(i-1);
end

table3 = [rh,r_ratio,eh,e_ratio];
T3 = array2table(table3)
%%
N = 63;
A = matPoisson(N);
f = fPossion(N);
Nmin = 3;
omega = 2/3;     
alpha1 = 2;
alpha2 = 1;
smoother = 1;
vold = zeros(N^2,1);

h = 1/(N+1);
x = 0:h:1;
y = 0:h:1;
temp = [];
for i = 1:N
    t = zeros(N,1);
    for j = 1:N
        t(j,1) = (x(i)^2-x(i)^4) * (y(j)^4-y(j)^2);
    end
    temp = [temp;t];
end

vexact = temp;

rh = zeros(15,1);
eh = zeros(15,1);

for i = 1:15
   vnew = mgVcycle(vold,f,alpha1,alpha2,omega,Nmin,N,smoother);
   e = vold - vnew;
   vold = vnew;
   r = A*e;
   eh(i) = norm(e);
   rh(i) = norm(r);
end

r_ratio = zeros(15,1);
e_ratio = zeros(15,1);

for i = 2:15
    r_ratio(i) = rh(i)/rh(i-1);
    e_ratio(i) = eh(i)/eh(i-1);
end

table4 = [rh,r_ratio,eh,e_ratio];
T4 = array2table(table4)
