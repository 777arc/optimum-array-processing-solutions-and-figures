%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 3.28															
% Beampattern of Dolph-Chebychev & Riblet								
% Xiaomin Lu  11/2/98
% K. Bell 9/5/00, Lillian Xu 04/16/2001, K. Bell 7/22/01, 9/30/01
%function used:  cheby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 21;										%21-element array
u = -1:1/200:1;

dB = 30;                                 %SLL= -30dB	
dRange = [1/2 1/4];				     	%d = 1/2, 1/4

m = 1;
for d = dRange
    R = 10^(dB/20);
    x0 = cosh( 2/(N-1)*acosh(R));
    alpha = 2*pi*d;
    x = 1/(1-cos(alpha))*( (x0+1)*cos(alpha*u)-1-x0*cos(alpha));
    beam = cheby((N-1)/2,x);
    beam = abs(beam)/max(abs(beam));
    table1(m,:) = 20*log10(beam);		%Riblet beampattern
    
    x0 = cosh( 1/(N-1)*acosh(R));
    x = x0*cos(pi*d*u);
    beam = cheby(N-1,x);
    beam = abs(beam)/max(abs(beam));
    table2(m,:) = 20*log10(beam);		%Chebychev beampattern
    
    m = m+1;
end


plot(u,table1(1,:),'-')
%text(-0.95,-5,'d=\lambda/2')
grid
hold on
plot(u,table2(1,:),'--');
h=legend('Riblet','Dolph-Chebychev');
axis([-1 1 -60 0])
%text(-0.025,-68,'(a)')
%title('Chebychev & Riblet beampattern, N=21, SLL=-20dB, d=\lambda/2')
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)

figure
plot(u,table1(2,:),'-')
%text(-0.95,-5,'d=\lambda/4')
grid
hold on
plot(u,table2(2,:),'--');
h=legend('Riblet','Dolph-Chebychev');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
axis([-1 1 -60 0])
%text(-0.025,-68,'(b)')
%title('Chebychev & Riblet beampattern, N=21, SLL=-20dB, d=\lambda/4')

