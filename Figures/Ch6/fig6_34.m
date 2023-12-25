%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.34
% K. Bell 7/26/01, 11/13/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

s = [0:0.01:1];   % sin^2(vm,va;rhon^-1)
M = [0.5 1 2 3 4 5 6 7 8 10 12 16];
nm = length(M);
ns = length(s);

A = zeros(nm,ns);

for m=1:nm
    A(m,:) = 1+(2*M(m)+M(m)^2)*s;
end

figure
plot(s,A);
xlabel('sin^{2}({\bfv}_{\itm}, {\bfv}_{\ita}; \rho_{\itn}^{-1})','Fontsize',14)
ylabel('{\itA}_{\itmvdr} /{\itA}_{\itmpdr}','Fontsize',14)
axis([0 1.4 0 80])

for m=1:nm
    I = find(A(m,:)<80);
    x = s(max(I));
    %y = A(m,max(I));
    y = min(A(m,ns),78);
    text(x+0.02,y-0.02,num2str(M(m)),'Fontsize',12)
end
