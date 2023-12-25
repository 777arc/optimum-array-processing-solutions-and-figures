%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.37
% layer noise power distribution
% K. Bell 7/24/01, 10/23/01
% function called: polarkb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

alpha = [1 0.5 0.25 0.125];
na = length(alpha);

theta = [0:0.01:1]*pi;
for n=1:na
    S = 1-alpha(n)/4-3*alpha(n)*cos(2*theta)/4;
    polarkb(theta,S)
    if n==1
        text(0.1,S(length(S)),['\alpha=1'],'Fontsize',12)
    else
        text(0.1,-S(length(S)),['\alpha=1/' int2str(1/alpha(n))],'Fontsize',12)
    end
    hold on
end
hold off
