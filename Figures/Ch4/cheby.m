%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cheby														
% computes Nth order chebychev polynomial T_N(x)								
% Xiaomin Lu  11/2/98
% Last updated 9/15/00 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T] = cheby(N,x)

T0 = 1;
T1 = x;
for n = 2 : N
   T = 2*x.*T1 - T0;
   T0 = T1;
   T1 = T;
end
