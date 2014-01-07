function [x] = modinv(a,b)
% x = modinv(a,b)
%
% Compute the inverse of a modulo b, given GCD(a,b) = 1
% using the Extended Euclidian Algorithm.

x1 = 1;
y1 = 0;
r1 = a;

x2 = 0;
y2 = 1;
r2 = b;

while (r2 ~= 0)
  r = mod(r1, r2);
  q = (r1 - r)/r2;

  x = x1 - q*x2;
  y = y1 - q*y2;

  r1 = r2;
  r2 = r;

  x1 = x2; 
  x2 = x;

  y1 = y2;
  y2 = y;
end

x = x1;
if (x < 0)
  x = x+b;
end

