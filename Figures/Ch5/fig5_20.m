%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.20
% Prolate spheroidal functions
% M-element array, d=lambda/2
% updated by Lillian Xu 04/05/2001
% K. Bell 7/23/01, 10/23/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
pack

M = 11;
d = 0.5;                        % sensor spacing wrt wavelength
x = sort([0.2*M/(2*pi) linspace(1e-2,2,200)]);   % deltau/BWNN
%D = [-(M-1)/2:1:(M-1)/2];   % sensor positions in wavelengths
D = [0:1:M-1];

% weights, normalized so that w(0)=1
p = [(2*pi/M)*x];      % (deltau/BWNN) -> p
np = length(p);
cutoff = [0.99 0.999 0.9999];
nco = length(cutoff);
siglam = zeros(nco,np);

for i = 1:nco
  cutoff(i);
  for n = 1:np
    A = toeplitz(sinc(D*p(n)))/M;
    [e,lam] = eig(A);
    [lam,ind] = sort(diag(abs(lam)));        % sort eigenvalues in ascending order 
    e = e(:,ind);                            % arrange eigenvectors in same order
    pwr = lam(M);
    m = M-1;
    while (m > 0) & (pwr <= cutoff(i))
      pwr = pwr+lam(m);
      m = m-1;
    end % while
    siglam(i,n) = M-m;
  end
end

%figure(1); stairs(x, siglam(1,:));
%figure(2); stairs(x, siglam(2,:));
%figure(3); stairs(x, siglam(3,:));

figure(1); clf;
for i = 1: nco
  subplot(nco,1,i);
  stairs(x, siglam(i,:));
  ylabel('{\itD}_{\its}({\itM})','Fontsize',14)
  title(['{\itM}=',num2str(10*log10(1/(1-cutoff(i)))),' dB'],'Fontsize',12);
end
  xlabel(['{\itu}_{\Delta}/{\itBW}_{\itNN}'],'Fontsize',14);
