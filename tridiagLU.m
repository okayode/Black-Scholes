function [e,f] = tridiagLU(a,b,c)

n = length(a);
e = zeros(n,1); f = e;
e(1) = b(1);
f(1) = c(1)/b(1);
for i=2:n
e(i) = b(i) - a(i)*f(i-1);
f(i) = c(i)/e(i);
end