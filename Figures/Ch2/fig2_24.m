%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.24
% HPBW versus steering angle
% Lillian Xu
% Modified by Xin Zhang
% Last updated by L. Xu 11/30/00, K. Bell 7/22/01, 10/4/01, 10/17/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

Nd_set=[1000:-5:700 700:-0.5:200 200:-0.05:20 15:-0.005:1];
theta_deg=[2.5 5 10 20 30 45 90];
%theta_deg=[2.5 5 10 20 30 45 60 90];
theta_set=theta_deg*pi/180;
figure
for num=1:size(theta_set,2)
    theta=theta_set(num);
    bw=0;
    for num1=1:size(Nd_set,2)
        Nd=Nd_set(num1);
        r=cos(theta)+0.443/Nd;
        if abs(r)<=1
            bw(num1)=acos(r)-acos(cos(theta)-0.443/Nd);
        end
    end
    loglog(Nd_set(1:size(bw,2)),abs(bw*180/pi));
    hold on
    rightpoint=abs(bw(1)*180/pi);
    text(1040,rightpoint,[num2str(theta_deg(num)),'^o'],'fontsize',10); % used to be 7
end
axis([1 1000 0.04 100])
xlabel('{\itNd}/\lambda','Fontsize',14)
ylabel('3-dB beamwidth in degrees','Fontsize',14)

scanlimit=-acos(1-2*0.443./Nd_set);
loglog(Nd_set,abs(scanlimit*180/pi),'--')

Endfire = 2*acos(1-0.433./Nd_set);
loglog(Nd_set,abs(Endfire*180/pi))
rightpoint=abs(Endfire(1)*180/pi);
text(1040,rightpoint,'\theta=0^o','fontsize',10); % used to be 7
%title('HPBW versus steering angle: standard linear array with uniform weighting')
text(12,35,'Endfire','Fontsize',12)
text(12,14.4,'Scan limit','Fontsize',12)
text(12,1.8,'Broadside','Fontsize',12)
plot([1 1000],[0.1 0.1],':')
plot([1 1000],[1 1],':')
plot([1 1000],[10 10],':')
plot([1 1000],[100 100],':')
plot([10 10],[0.04 200],':')
plot([100 100],[0.04 200],':')
axis([1 1000 0.04 200])
