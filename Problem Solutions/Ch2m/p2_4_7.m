%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.4.7
% Xin Zhang
% Last updated 2/11/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

l = 1;
for N = 2:30

n = (-(N-1)/2:(N-1)/2).';
u = (-1:0.001:1);
vv = exp(j*n*pi*u);
w = 1/N*[ones(N,1)];
B = w'*vv;
B = abs(B).^2;

%Taylor Approximation
if N == 10 | N == 11 | N == 21
   B2 = 1-1/12*(N^2-1)*pi*pi*(u.^2);
   if N == 10
      figure;
      plot(u,B,'-',u,B2,'--')
 	   grid on
 	   ylabel('B^2')
      xlabel('u')
      legend('exact beampattern','approximation')
      title(['Problem 2.4.7, N = ',num2str(N)])
      axis([0 0.2 0 1])
      figure;
   end
   subplot(3,1,l);
   plot(u,B,'-',u,B2,'--')
   l = l+1;
   grid on
   ylabel('B^2')
   if N == 21
      xlabel('u')
   end
   title(['N = ',num2str(N)])
   if N == 10
      title(['Problem 2.4.7, N = ',num2str(N)])
      legend('exact beampattern','approximation')
   end
   axis([-1 1 0 1])
end

% Find HPBW
u = [0:0.0001:0.2];
vv = exp(j*n*pi*u);
Bs = w'*vv;
Bs = abs(Bs).^2-0.5;
Bs2 = 0.5-1/12*(N^2-1)*pi*pi*(u.^2);
[y,I] = min(abs(Bs));
umin = u(I);
HPBW_u = 2*umin;
x(N-1) = umin*N*pi/2;
Bmin = y;
[y,I] = min(abs(Bs2));
umin = u(I);
HPBW_u = 2*umin;
x2(N-1) = umin*N*pi/2;
Bmin = y;

end
N = (2:1:30);
%N,exact N*pi*d*u/lambda,Taylor approximation 
[N.' x.' x2.']