%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bpsphcut_theta.m 
% plots beampattern cuts for specified values of theta
% K. Bell 9/6/00
% updated by K. Bell 10/10/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bpsphcut_theta(D,w,lim,theta)

% D = Nx3 matrix of sensor locations
% w = Nx1 vector of weights
% lim = dB lower limit
% theta  = 2x1 vector of theta cut values

% reasonable values
% lim = -60;
% theta_deg = [90 60 30];
% theta = 2*pi*theta_deg/360;

phi = pi*[-1:0.001:1];
lp = length(phi);
lt = length(theta);

xx = zeros(lp,lt);
yy = zeros(lp,lt);
zz = zeros(lp,lt);
figure
for nth = 1:1:lt
   B = zeros(1,lp);
  for nph = 1:1:lp
     ux = sin(theta(nth))*cos(phi(nph));
     uy = sin(theta(nth))*sin(phi(nph));
     uz = cos(theta(nth));

     V = exp(j*2*pi*D*[ux;uy;uz]);
     B(nph) = 10*log10(abs(w'*V).^2);
  end
  subplot(lt,1,nth)
  plot(phi/pi,B)
  grid on
  xlabel('\phi/\pi')
  ylabel('Beampattern (dB)')
       title(['\theta = ' int2str(theta(nth)*180/pi) ' deg.'])

  axis([-1 1 -60 0])

end
