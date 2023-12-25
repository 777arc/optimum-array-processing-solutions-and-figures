%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.4, Example 7.3.1
%  K. Bell 1/16/03
% Updated by K. Bell 1/30/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N = 10;
d = 0.5;
L = 200;
K = [10 20 30 40 60 100 200 500 1000];
nk = length(K);
D  = [-(N-1)/2:1:(N-1)/2].';

us = 0;
sigma_s = 10^(10/10);
uI1 = 0.15;
sigma_1 =10^(10/10);

sigma_n = 1;

vs  = exp(j*2*pi*d*D*us);
vI1  = exp(j*2*pi*d*D*uI1);


Rs = sigma_s*vs*vs';
Rn =  sigma_1*vI1*vI1'+ sigma_n*eye(N);
Rninv = inv(Rn);

Wq = vs/N;            % C*inv(C'*C)*f  

Wo = Rninv*vs*inv(vs'*Rninv*vs);
SINRopt = real(Wo'*Rs*Wo)/real(Wo'*Rn*Wo);      


Erho = zeros(1,nk);
rho_bar = zeros(1,nk);
Erhofb = zeros(1,nk);
rho_barfb = zeros(1,nk);
for k = 1:nk
   rho = zeros(1,L);
   rhofb = zeros(1,L);
   for ll = 1:L
      ii1 = sqrt(sigma_1/2)*randn(1,K(k))+j*sqrt(sigma_1/2)*randn(1,K(k));  % interference samples
      nn  = sqrt(sigma_n/2)*randn(N,K(k))+j*sqrt(sigma_n/2)*randn(N,K(k));    % noise samples
      nv =  vI1*ii1+nn;
      nvb = conj(flipud(nv));
      
      Rnhat   = nv*nv'/K(k);
      Rnhatfb = 0.5*Rnhat + 0.5*(nvb*nvb')/K(k);
            
      Rnhatinv   = inv(Rnhat);
      Rnhatfbinv = inv(Rnhatfb);
     
      wv   = Rnhatinv*vs*inv(vs'*Rnhatinv*vs);
      wvfb = Rnhatfbinv*vs*inv(vs'*Rnhatfbinv*vs);
      rho(ll)   = (real(wv'*Rs*wv)/real(wv'*Rn*wv))/SINRopt;
      rhofb(ll) = (real(wvfb'*Rs*wvfb)/real(wvfb'*Rn*wvfb))/SINRopt;
   end % trials
   rho_bar(k) = sum(rho)/L;
   Erho(k) = (K(k)+2-N)/(K(k)+1);
   rho_barfb(k) = sum(rhofb)/L;
   Erhofb(k) = (2*K(k)+2-N)/(2*K(k)+1);
end

figure(1)
h4=semilogx(K,10*log10(Erho),'*');
hold on
h3=semilogx(K,10*log10(rho_bar),'-');
%h2=semilogx(K,10*log10(Erhofb),'db');
h1=semilogx(K,10*log10(rho_barfb),'-.');
hold off
legend([h3 h1 h4],'MVDR','MVDR.FB','{\it E}[\rho]',4)
grid on
xlabel('Number of snapshots (K)')
ylabel('Normalized average {\it SINR_0} (dB)')
