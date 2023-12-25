%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 5.3.1
% K. Bell 10/23/98
% updated by K. Bell 11/20/03
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dp = [0:0.01:5];
SR = sinc(2*dp);

figure(1)
clf
subplot(3,1,1)
plot(dp,SR)
%xlabel('|\Delta p|/\lambda')
ylabel('Re[S_x(\omega_o,\Delta_p)]')
title('Problem 5.3.1')
grid on
text(0.75,0.5,'Re[S_x(\omega_o,\Delta_p)] does not depend on \alpha or \theta_p')


alpha = 0.5;
theta_p = [0 pi/6 pi/3 pi/2].';
SI = alpha*cos(theta_p)*((sinc(2*dp)-cos(2*pi*dp))./(2*pi*dp));
%SI = ones(5,1)*sinc(2*dp);
%SI = ones(5,1)*cos(2*pi*dp);
subplot(3,1,2)
h1=plot(dp,SI(1,:),'-');
hold on
h2=plot(dp,SI(2,:),':');
h3=plot(dp,SI(3,:),'-.');
h4=plot(dp,SI(4,:),'--');
hold off
legend([h1 h2 h3 h4],['\theta_p = ' num2str(180*theta_p(1)/pi) ' deg.'],...
   ['\theta_p = ' num2str(180*theta_p(2)/pi) ' deg.'],...
   ['\theta_p = ' num2str(180*theta_p(3)/pi) ' deg.'],...
   ['\theta_p = ' num2str(180*theta_p(4)/pi) ' deg.']);
%xlabel('|\Delta p|/\lambda')
ylabel('Im[S_x(\omega_o,\Delta_p)]')
title(['\alpha = ' num2str(alpha)])
axis([min(dp) max(dp) -0.1 0.5])
text(0.75,0.27,'Plots change sign for \theta_p > 90 deg.')
grid on

alpha = [0 0.25 0.5 1].';
theta_p = pi/3;
SI = alpha*cos(theta_p)*((sinc(2*dp)-cos(2*pi*dp))./(2*pi*dp));
subplot(3,1,3)
h1=plot(dp,SI(1,:),'-');
hold on
h2=plot(dp,SI(2,:),':');
h3=plot(dp,SI(3,:),'-.');
h4=plot(dp,SI(4,:),'--');
hold off
legend([h1 h2 h3 h4],['\alpha = ' num2str(alpha(1))],['\alpha = ' num2str(alpha(2))],...
   ['\alpha = ' num2str(alpha(3))],['\alpha = ' num2str(alpha(4))]);
xlabel('|\Delta p|/\lambda')
ylabel('Im[S_x(\omega_o,\Delta_p)]')
title(['\theta_p = ' num2str(180*theta_p/pi) ' deg.'])
grid on
  axis([min(dp) max(dp) -0.1 0.5])
 
set(gcf,'Paperposition',[0.25 1 8 9])


