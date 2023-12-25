%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bpsphere.m plots 3D beampattern
% K. Bell
% 9/6/00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bpsphere(D,w,lim,m)

% D = Nx3 matrix of sensor locations
% w = Nx1 vector of weights
% lim = dB lower limit
% m = points evaluated
% some reasonable values
%lim = -60;
%m=100;

% -pi <= phi <= pi is a row vector.
% 0 <= theta <= pi is a column vector.

phi = (-m:2:m)/m*pi;
theta = (0:1:m)'/m*pi;
lp = length(phi);
lt = length(theta);

xx = zeros(lp,lt);
yy = zeros(lp,lt);
zz = zeros(lp,lt);

for nph = 1:1:lp
  for nth = 1:1:lt
     ux = sin(theta(nth))*cos(phi(nph));
     uy = sin(theta(nth))*sin(phi(nph));
     uz = cos(theta(nth));

     V = exp(j*2*pi*D*[ux;uy;uz]);
     B = max(lim,10*log10(abs(w'*V).^2));
     B = -lim+B;
     xx(nph,nth) = B*ux;
     yy(nph,nth) = B*uy;
     zz(nph,nth) = B*uz;
  end
end
figure
clf
surf(xx,yy,zz);
axis((-lim+10)*[-1 1 -1 1 -1 1])
view([-30 10])
hold on
h=plot3([1 -1]*lim-10,[0 0],[0 0],'k');
set(h,'Linewidth',2)
h=plot3([0 0],[0 0], [1 -1]*lim-5,'k');
set(h,'Linewidth',2)
h=plot3([0 0], [1 -1]*lim-5,[0 0],'k');
set(h,'Linewidth',2)
hold off    
xlabel('x')
ylabel('y')
zlabel('z')

