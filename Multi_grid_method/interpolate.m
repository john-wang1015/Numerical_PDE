function I2hh = interpolate(N,N2h)
% Builds the iterpolation marix from level 2h to level h.
% There are N2h*N2h interior points on the coarse level (level 2h), 
% and N*N interior points on the fine level (level h).
% This means that I2hh is a (N*N) * (N2h*N2h) matrix (many more
% rows than columns). I2hh is again a sparse matrix
I2hh = spalloc(N^2,N2h^2,10*N^2*N^2);
steretl = [1/4 1/2 1/4 1/2 1 1/2 1/4 1/2 1/4];

for i = 1:N2h
    for j = 1:N2h
        k2h = (j-1)*N2h+i;
        sw = (2*j-2)*N + 2*i-1; %
        s = (2*j-2)*N + 2*i;
        se = (2*j-2)*N + 2*i+1;
        w = (2*j-1)*N + 2*i-1;
        c = (2*j-1)*N + 2*i;
        e = (2*j-1)*N + 2*i +1;%
        nw = (2*j)*N + 2*i-1;
        n = (2*j)*N + 2*i;
        ne = (2*j)*N + 2*i+1;
        index = [sw,s,se,w,c,e,nw,n,ne];
        I2hh(index,k2h) = steretl;
    end
end

end
