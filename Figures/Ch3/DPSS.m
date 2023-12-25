%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DPSS.m
% calculates DPSS weights
% Lillian Xiaolan Xu
% Last updated 9/6/00 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w=DPSS(psi0,N)
	for ma=1:N
		for na=1:N
			A(ma,na)=2*psi0*sinc((ma-na)*psi0/pi);		
		end;
	end;

	[E eigv]=eig(A);
	[Y,index]=sort(real(diag(eigv)));
	w=E(:,index(N));
				%Discrete Prolate Spheroidal Function
