%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% taylor_circ														
%   Taylor distribution < g(p) > for circular aperture
% Xiaomin Lu  11/2/98
% Last updated 10/11/00 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function weight = taylor_circ(p)

Bessel_null =    [1.21966989126650;2.23313059438153;3.23831548416624;...
   4.24106286379607;5.24276437687019;6.24392168986449];
nBar = 6;
SLL = 20;
scale = 5;

A = 1/pi*acosh(10^(SLL/20));
weight = 1 ;

for m = 1:nBar-1
      F = -besselj(0,pi*Bessel_null(m));
      for n = 1:nBar-1
         
         Zn = Bessel_null(nBar)*sqrt((A^2+(n-0.5)^2)./(A^2+(nBar-0.5)^2));
         F = F*(1-Bessel_null(m)^2/Zn^2);
         if (n ~= m) F = F/(1-Bessel_null(m)^2/Bessel_null(n)^2); end
      end
      weight = weight+F/besselj(0,pi*Bessel_null(m))^2*besselj(0,Bessel_null(m)*p);

end
weight = weight*2/pi^2;

   




