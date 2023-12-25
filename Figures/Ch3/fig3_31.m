%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.31
% Beam pattern of Villeneuve weighting
% 	 (a) N = 21
% 	 (b) N = 41
% Lillian Xu
% Modified by Xin Zhang 3/24/99, K. Bell 9/5/00
% Lillian Xu 04/16/2001, K. Bell 7/22/01, 9/30/01
% function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
N=21;
M=(N-1)/2;
n1=6;
d=1/2;

sidelobe=20;    % *dB below the main lobe maximum
R=10.^(sidelobe/20);
x0=cosh(1/2/M*acosh(R));
u=-1:0.001:1;

%--------------------
c=pi*n1/N/acos(1/x0*cos((2*n1-1)*pi/4/M));

b2=sinc(N/2*u)./sinc(u/2);
for num=1:n1-1
	up=2*acos(1/x0*cos((2*num-1)*pi/4/M))/pi;
	b2=b2.*sin(pi/2*(u-c*up)).*sin(pi/2*(u+c*up))./sin(pi/2*(u-2*num/N))./sin(pi/2*(u+2*num/N));
end
for num=1:size(b2,2)
	if b2(num)==Inf | b2(num)==-Inf
		if num~=1 & num~=size(b2,2)
			b2(num)=(b2(num+1)+b2(num-1))/2;
		elseif num==1
			b2(num)=b2(num+1);
		else
			b2(num)=b2(num-1);
		end
	end
end

b2=b2./max(b2);

figure
plot(u,20*log10(abs(b2)));
axis([0 1 -50 0]);
%title(['Beampatterns of Villeneuve Distribution, N=',num2str(N), ...
%      ', sidelobe=',num2str(sidelobe),'dB, n1=',num2str(n1) ])
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
grid on
hold off

%**** Figure (b)
N=41;
M=(N-1)/2;
x0=cosh(1/2/M*acosh(R));

%--------------------
c=pi*n1/N/acos(1/x0*cos((2*n1-1)*pi/4/M));

b2=sinc(N/2*u)./sinc(u/2);
for num=1:n1-1
	up=2*acos(1/x0*cos((2*num-1)*pi/4/M))/pi;
	b2=b2.*sin(pi/2*(u-c*up)).*sin(pi/2*(u+c*up))./sin(pi/2*(u-2*num/N))./sin(pi/2*(u+2*num/N));
end
for num=1:size(b2,2)
	if b2(num)==Inf | b2(num)==-Inf
		if num~=1 & num~=size(b2,2)
			b2(num)=(b2(num+1)+b2(num-1))/2;
		elseif num==1
			b2(num)=b2(num+1);
		else
			b2(num)=b2(num-1);
		end
	end
end

b2=b2./max(b2);

figure
plot(u,20*log10(abs(b2)));
axis([0 1 -50 0]);
%title(['Beampatterns of Villeneuve Distribution, N=',num2str(N), ...
%      ', sidelobe=',num2str(sidelobe),'dB, n1=',num2str(n1) ])
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
grid on
hold off