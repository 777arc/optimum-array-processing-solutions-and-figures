%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 3.4.15
% K. Bell 10/6/99
% Last Updated 9/02/02 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%Na=[13 21 51 75 101 151 201 251 301 351 401 451 501];
Na=[13 21 31 51 63 69 75 101 201 301 351 401 501 701 1000];
Na=[10 20 30 50 70 75 100 200 300 350 400 500 700 800 1000];
%Na=[13 21 31 51 63 65];
nN = length(Na);
SL = [-15 -20 -30 -40 -60];
l_style = ['-d';'-+';'-*';'-o';'-s'];
nS = length(SL);
DID = zeros(nS,nN);

d=0.5;


for q=1:nN
   N = Na(q);
   for k=1:nS
      R = 10^(-SL(k)/20);
      x0 = cosh(acosh(R)/(N-1));
      uk = [-(N-1)/2:1:(N-1)/2]*1/(N*d);
      x = x0*cos(pi*d*uk);
      % Alternate expression for Tm(x) from Gradshteyn & Ryzhik,8.940,p.1032
      s= j*sqrt(1-x.^2);
      Bp = 0.5*((x+s).^(N-1) + (x-s).^(N-1));
      Bp = Bp/R;
      Bk = conj(Bp).*exp(-j*2*pi*d*uk*(N-1)/2);
      b = ifft(Bk);
      w = b.*exp(-j*[0:1:N-1]*pi*((N-1)/2)*2/N);     % w are complex with imag. part = 0
      w = real(w.');                                 % gets rid of zero imag. part
      w=w/sum(w);                                    % normalization not really necessary
      DID(k,q) = -10*log10(w'*w);
   end
end

figure(1)
clf
for k= 1:5
   h(k) = semilogx(Na,DID(k,:),l_style(k,:));
   hold on
end
hold off
   xlabel('N')
   ylabel('Directivity Index')
   title(['Problem 3.4.15, Dolph-Chebychev'])
   set(gca,'XTick',[10 100 1000])
   set(gca,'XTickLabel',['  10';' 100';'1000'])
   legend(h,'-15 dB','-20 dB','-30 dB','-40 dB','-60 db',4)
set(gcf,'Paperposition',[0.25 4.5 8 6])


