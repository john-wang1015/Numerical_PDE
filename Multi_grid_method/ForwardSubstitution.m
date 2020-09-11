function x = ForwardSubstitution(Q, b)
x = b;

for r = 2:length(b)
  for c = 1:r-1
    x(r) = x(r) - x(c) * Q(r,c);
  end
x = Q\b;
end
