%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.3
%  K. Bell 1/22/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N = 10;
d = 0.5;
K = [1;2;6;10;15]*N;
nk = length(K);
rho = [0:0.01:1];
figure(1)
for k=1:nk
    %gamma(n+1)=n!
    
    p = (gamma(K(k)+1)/(gamma(N-1)*gamma(K(k)+2-N)))*((1-rho).^(N-2)).*(rho.^(K(k)+1-N));
    switch k
        case 1 
            ls = '-';
        case 2 
            ls = '-.';
        case 3 
            ls = '--';
        case 4 
            ls = ':';
        case 5 
            ls = ':x';
    end
    h4(k)=plot(rho,p,ls);
    hold on
end

legend([h4],'{\it K}={\it N}','{\it K}=2{\it N}','{\it K}=6{\it N}','{\it K}=10{\it N}','{\it K}=15{\it N}',2)
xlabel('\rho')
ylabel('Probability density')
hold off