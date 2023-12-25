%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sinc.m
% same as Matlab's sinc function
% sinc(x) =  sin(pi*x)/(pi*x)
% Last updated 8/31/00 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=sinc(x)
y=ones(size(x));
i=find(x);
y(i)=sin(pi*x(i))./(pi*x(i));
