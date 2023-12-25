function w = chebklb(N,SL)

R = 10^(SL/20);
x0 = cosh(acosh(R)/(N-1));
p = [1:1:N-1];
a = 2*acos(cos((2*p-1)*pi/(2*(N-1)))/x0);
z = exp(j*a);
wt = real(poly(z).');

w = wt/max(wt);
