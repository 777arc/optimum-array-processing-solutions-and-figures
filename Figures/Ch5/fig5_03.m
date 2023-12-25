%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.3
% Functions in integrand 
% (a) sinc^2() functions for various m
% (b) spectrum of bandpass process (B_s Delta_T = 16)
% Xin Zhang 6/20/00
% Updated by K. Bell 10/12/00, 10/23/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

w_wd = -8:0.01:8;
m_set = -8:1:8;
figure
for l = 1:length(m_set)
   m = m_set(l);
   y1(l,:) = sinc(w_wd-m).^2;
   plot(w_wd,y1(l,:));
   hold on
end
text(-8,1.1,'{\itm}=-8','FontSize',12,'HorizontalAlignment','Left')
text(8,1.1,'{\itm}=8','FontSize',12,'HorizontalAlignment','Right')
for m=-6:2:6
    text(m,1.1,['{\itm}=' int2str(m)],'HorizontalAlignment','Center','Fontsize',12)
end
for m=-7:2:7
    text(m,1.03,['{\itm}=' int2str(m)],'HorizontalAlignment','Center','Fontsize',12)
end


hold off
xlabel('{\it\omega}_{\itL}/{\it\omega}_{\Delta}','FontSize',14)
axis([-8 8 0 1.3])

B_d = 16;
w = -10:0.01:10;
height = 1;
for l = 1:length(w)
   if (w(l)>=-B_d/2)&(w(l)<=B_d/2)
      spec(l)=height;
   else
      spec(l)=0;
   end
end
figure
plot(w,spec)
xlabel('{\it\omega}_{\itL}/{\it\omega}_{\Delta}','FontSize',14)
%ylabel('[{S}_{\itx}({\omega}_{L}+{\omega}_{c})]_{nn}/(\sigma_{s}^2/{B}_{s})','FontSize',14)
ylabel('[{\itS}_{x}({\it\omega}_{\itL}+{\it\omega}_{\itc})]_{\itnn}/(\it\sigma_{\its}^2/{\itB}_{\its})','FontSize',14)
title('{\itB}_{\its} {\it\DeltaT} = 16','Fontsize',12)
axis([-10 10 0 5])
