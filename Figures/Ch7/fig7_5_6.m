%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.5, 7.6 Example 7.3.2
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
SNR = [0:10:30];
ns = length(SNR);
uI1 = 0.29;
uI2 = 0.45;
sigma_1 =10^(20/10);
sigma_2 =10^(20/10);

sigma_n = 1;

vs  = exp(j*2*pi*d*D*us);
vI1  = exp(j*2*pi*d*D*uI1);
vI2  = exp(j*2*pi*d*D*uI2);

Rn =  sigma_1*vI1*vI1'+sigma_2*vI2*vI2'+ sigma_n*eye(N);
Rninv = inv(Rn);
Wq = vs/N;            % C*inv(C'*C)*f  
Wo = Rninv*vs*inv(vs'*Rninv*vs);


rho_bar = zeros(ns,nk);
rho_barfb = zeros(ns,nk);
eta_bar = zeros(ns,nk);
eta_barfb = zeros(ns,nk);

for n = 1:ns
    sigma_s = 10^(SNR(n)/10);
    Rs = sigma_s*vs*vs';
    Rx = Rs+Rn;
    
    SINRopt = real(Wo'*Rs*Wo)/real(Wo'*Rn*Wo);      
        
    Erho = zeros(1,nk);
    Eeta1 = zeros(1,nk);
    for k = 1:nk
        rho = zeros(1,L);
        rhofb = zeros(1,L);
        eta = zeros(1,L);
        etafb = zeros(1,L);
        for ll = 1:L
            ss  = sqrt(sigma_s/2)*randn(1,K(k))+j*sqrt(sigma_s/2)*randn(1,K(k));  % signal samples
            ii1 = sqrt(sigma_1/2)*randn(1,K(k))+j*sqrt(sigma_1/2)*randn(1,K(k));  % interference samples
            ii2 = sqrt(sigma_2/2)*randn(1,K(k))+j*sqrt(sigma_2/2)*randn(1,K(k));  % interference samples
            nn  = sqrt(sigma_n/2)*randn(N,K(k))+j*sqrt(sigma_n/2)*randn(N,K(k));    % noise samples
            
            sv =  vs*ss;
            nv =  vI1*ii1+vI2*ii2+nn;
            xx  = sv+nv;
            
            nvb = conj(flipud(nv));
            xxb = conj(flipud(xx));
            
            Rnhat   = nv*nv'/K(k);
            Rnhatfb = 0.5*Rnhat + 0.5*(nvb*nvb')/K(k);
            Rxhat   = xx*xx'/K(k);
            Rxhatfb = 0.5*Rxhat + 0.5*(xxb*xxb')/K(k);
            
            Rnhatinv   = inv(Rnhat);
            Rnhatfbinv = inv(Rnhatfb);
            Rxhatinv   = inv(Rxhat);
            Rxhatfbinv = inv(Rxhatfb);
            
            wv   = Rnhatinv*vs*inv(vs'*Rnhatinv*vs);
            wvfb = Rnhatfbinv*vs*inv(vs'*Rnhatfbinv*vs);
            wp   = Rxhatinv*vs*inv(vs'*Rxhatinv*vs);
            wpfb = Rxhatfbinv*vs*inv(vs'*Rxhatfbinv*vs);

            rho(ll)   = (real(wv'*Rs*wv)/real(wv'*Rn*wv))/SINRopt;
            rhofb(ll) = (real(wvfb'*Rs*wvfb)/real(wvfb'*Rn*wvfb))/SINRopt;
            eta(ll)   = (real(wp'*Rs*wp)/real(wp'*Rn*wp))/SINRopt;
            etafb(ll) = (real(wpfb'*Rs*wpfb)/real(wpfb'*Rn*wpfb))/SINRopt;
        end % trials
        rho_bar(n,k) = sum(rho)/L;
        Erho(k) = (K(k)+2-N)/(K(k)+1);
        rho_barfb(n,k) = sum(rhofb)/L;
        eta_bar(n,k) = sum(eta)/L;
        a=K(k)-N+2;
        b = N-1;
        Eeta1(k) = (a/(a+b))/(1+SINRopt*b/(a+b+1));
        eta_barfb(n,k) = sum(etafb)/L;
    end
    
    figure
    h1=semilogx(K,10*log10(rho_bar(n,:)*SINRopt),'-');
    hold on
    h2=semilogx(K,10*log10(rho_barfb(n,:)*SINRopt),'-.');
    h3=semilogx(K,10*log10(eta_bar(n,:)*SINRopt),'--');
    h4=semilogx(K,10*log10(eta_barfb(n,:)*SINRopt),':');
    h5=semilogx(K,10*log10(Erho*SINRopt),'o');
    h6=semilogx(K,10*log10(Eeta1*SINRopt),'x');
    hold off
    legend([h1 h2 h3 h4 h5 h6],'MVDR','MVDR-FB','MPDR','MPDR-FB','MVDR-analytic','MPDR-analytic')
    grid on
    xlabel('Number of snapshots (K)')
    ylabel('Average {\it SINR_0} (dB)')
    axis([10 1000 -10 40])
end

figure
semilogx(K,10*log10(eta_bar(1,:)),'-');
hold on
semilogx(K,10*log10(eta_bar(2,:)),'-.');
semilogx(K,10*log10(eta_bar(3,:)),'--');
semilogx(K,10*log10(eta_bar(4,:)),':');
semilogx(K,10*log10(eta_barfb(1,:)),'-x');
semilogx(K,10*log10(eta_barfb(2,:)),'-.o');
semilogx(K,10*log10(eta_barfb(3,:)),'--+');
semilogx(K,10*log10(eta_barfb(4,:)),':*');
hold off
legend('{\it SNR} = 0 dB','{\it SNR} = 10 dB','{\it SNR} = 20 dB','{\it SNR} = 30 dB','{\it SNR} = 0 dB, FB','{\it SNR} = 10 dB, FB','{\it SNR} = 20 dB, FB','{\it SNR} = 30 dB, FB')
grid on
xlabel('Number of snapshots (K)')
ylabel('Normalized average {\it SINR_0} (dB)')
axis([10 1000 -50 0])
