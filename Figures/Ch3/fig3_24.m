%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.24
% Dolph-Chebychev Polynomial
% Lillian Xu 
% Zhang Xin 3/23/99, Lillian Xu 04/16/2001, K. Bell 7/22/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

x = -1.2:0.0001:1.2;
T5 = 16*x.^5-20*x.^3+5.*x;
x0 = 1.2;
R = 16*x0^5-20*x0^3+5*x0;

l = 0;
for k = 1:length(x)
   if ((x(k) > 0.8) & (x(k) < 1) & (abs(T5(k)) < 0.0008))
      l = x(k);				% l = 0.9511;
   end
end

figure
plot(x,T5)
hold on;
plot(x0*[-1 1],[1 1],'--',x0*[-1 1],[-1 -1],'--');
plot(x0*[1 1],[-12.3 12.5],'--');
plot(x0,R,'.');
plot(l*[1,1],[-5.4,0],'--');
plot((-1)*l*[1,1],[-5.4,0],'--');
hold off;
text(x0-0.03,-13,'{\it x}_{0}','Fontsize',12);
text(x0+0.02,R,'{\it( x}_{0},{\it R)}','Fontsize',12);
text(l-0.02,-5.8,'{\it x}_{1}','Fontsize',12);
text(-l-0.02,-5.8,'{\it -x}_{1}','Fontsize',12);
grid
xlabel('\it x','Fontsize',14)
ylabel('{\it T}_5{\it(x)}','Fontsize',14)
%title('Dolph-Chebychev Polynomial')
