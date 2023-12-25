%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 3.7.5
% K. Bell 10/19/99
% Updated 9/28/00 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=41;
nbar = 6;
SL = -30;
m=N-1;
d=0.5;
D = d*[-(N-1)/2:1:(N-1)/2].';
u = [-1:0.001:1];
v = exp(j*2*pi*D*u);
nu = length(u);

wn = zeros(N,3);


% villeneuve weights

R = 10^(-SL/20);
x0 = cosh(acosh(R)/m);
cc=cos((2*[1:1:floor((N-1)/2)]-1)*0.5*pi/(N-1)).';   % Dolph-chebychev roots
udc = acos(cc/x0)/(pi*d);

if rem(N,2)==0                   % N even
   uu = [1:1:(N/2)-1].'/(N*d)    % uniform roots
   u_nbar = udc(nbar);       % u_nbar Chebychev
   sigma = nbar/(N*d*u_nbar);
   udcmod = udc*sigma;           % modified Chebychev
   uz = [udcmod(1:nbar-1);uu(nbar:(N/2)-1)];
   ut = [uz;-uz;-1];
else                             % Nodd
   uu = [1:1:(N-1)/2].'/(N*d);   % uniform roots
   u_nbar = udc(nbar) ;      % u_nbar Chebychev
   sigma = nbar/(N*d*u_nbar);
   udcmod = udc*sigma;           % modified Chebychev
   uz = [udcmod(1:nbar-1);uu(nbar:(N-1)/2)];
   ut = [uz;-uz];
end

wv = poly(exp(j*2*pi*d*ut)).';
wv = wv/sum(wv);


% 0 order
%order = 0;
%null1 = [0.24 0.3 0.36];
%C1 = exp(j*2*pi*D*null1);
%PC1 = eye(N) - C1*inv(C1'*C1)*C1';

%null2 = [0.24 0.27 0.3 0.33 0.36];
%C2 = exp(j*2*pi*D*null2);
%PC2 = eye(N) - C2*inv(C2'*C2)*C2';

%null3 =[0.24 0.26 0.28 0.3 0.32 0.34 0.36];
%C3 = exp(j*2*pi*D*null3);
%PC3 = eye(N) - C3*inv(C3'*C3)*C3';

% 0 and 1 order
%order = 1;
%null1 = 0.3;
%C1 = [exp(j*2*pi*D*null1)  j*2*pi*D.*exp(j*2*pi*D*null1)] ;
%PC1 = eye(N) - C1*inv(C1'*C1)*C1';

%null2 = [0.27 0.33];
%C2 = [exp(j*2*pi*D*null2)  j*2*pi*D*[1 1].*exp(j*2*pi*D*null2)] ;
%PC2 = eye(N) - C2*inv(C2'*C2)*C2';

%null3 = [0.26 0.3 0.34];
%C3 = [exp(j*2*pi*D*null3)  j*2*pi*D*[1 1 1].*exp(j*2*pi*D*null3)] ;
%PC3 = eye(N) - C3*inv(C3'*C3)*C3';

% 0, 1, and 2 order
order = 2;
null1 = 0.3;
C1 = [exp(j*2*pi*D*null1)  j*2*pi*D.*exp(j*2*pi*D*null1) ((j*2*pi*D).^2).*exp(j*2*pi*D*null1)] ;
PC1 = eye(N) - C1*inv(C1'*C1)*C1';

null2 = [0.27 0.33];
C2 = [exp(j*2*pi*D*null2)  j*2*pi*D*[1 1].*exp(j*2*pi*D*null2)...
   ((j*2*pi*D*[1 1]).^2).*exp(j*2*pi*D*null2)] ;
PC2 = eye(N) - C2*inv(C2'*C2)*C2';

null3 = [0.26 0.3 0.34];
C3 = [exp(j*2*pi*D*null3)  j*2*pi*D*[1 1 1].*exp(j*2*pi*D*null3) ...
((j*2*pi*D*[1 1 1]).^2).*exp(j*2*pi*D*null3)] ;
PC3 = eye(N) - C3*inv(C3'*C3)*C3';

w(:,1) = PC1*wv;
w(:,1) = w(:,1)/sum(w(:,1));

w(:,2) = PC2*wv;
w(:,2) = w(:,2)/sum(w(:,2));

w(:,3) = PC3*wv;
w(:,3) = w(:,3)/sum(w(:,3));

for qq=1:3
   Bw = w(:,qq)'*v;
   figure(1)
   subplot(3,1,qq)
   plot(u,10*log10(abs(Bw).^2),'r')
   hold on
   plot([0.24 0.24],[-60 0],'--')
   plot([0.36 0.36],[-60 0],'--')
   hold off
   xlabel('u=\psi/\pi')
   ylabel('Beampattern (dB)')
   set(gca,'YTick',[-60:10:0])   
   if qq ==1
      title(['Problem 3.7.5, ' int2str(length(null1)) ' ' int2str(order) '-order nulls at u= ' num2str(null1)])
   elseif qq==2
      title([int2str(length(null2)) ' ' int2str(order) '-order nulls at u= ' num2str(null2)])
   else
      title([int2str(length(null3)) ' ' int2str(order) '-order nulls at u= ' num2str(null3)])
   end
   axis([-1 1 -60 0])
   
end
set(gcf,'Paperposition',[0.25 1 8 9])