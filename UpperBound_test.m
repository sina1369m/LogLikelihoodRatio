clear;
clc;
close all


r=1;
theta=[1];
theta=0.6*theta/(sum(theta));

%Generate Data
n=10^5;
%Find the upper bound
%Param
ccc=0;
cc2=0
for c2=0.001:0.001:0.5  
    cc2=cc2+1;
    cc1=0
    for c1=-1:0.001:-0.001    
        cc1=cc1+1;
        ccc=ccc+1;
        delp=n^(c1);
        del=n^(c2);
        tr=1-sum(theta);%Theta_res
        tm=min(min(theta),tr); %Theta_min 
        q=max((r+1)/tm, r*(tm^(-2)-tm^(-1))-1);
        w=max(2*r, tm^(-2)-tm^(-1) +r*(4+r*( tm^(-2)-tm^(-1) )))/tm;
        %h=400*r^0.25*(1/tm -1)^1.5; 
        h=400*(r^0.25)*(sum(theta.*(1./theta -1).^1.5) + tr*(1/tr -1)^1.5); 
        eta=r^5/(tm^6);
        zeta=r^3/(3*tm^3);
        dels=delp/(1-delp);
        a=1;%CHOOSE the worse a later
        Delp=a*dels + (eta*del^2)/(n*(1-delp))+(zeta*del^1.5)/sqrt(n);

        eps(ccc)=2*r*exp(-0.5*n*(delp^2)/(w+delp*q/3));
        eps2(ccc)=2*r*exp(-2*del/r);
        eps3(ccc)=h/sqrt(n);
        eps4(ccc)=Delp*(a/2 + Delp)^(0.5*r-1);
        UB(cc2,cc1)=eps(ccc)+eps2(ccc)+eps3(ccc)+eps4(ccc);
    end
end
C2=repmat(0.001:0.001:0.5,cc1,1);
C1=repmat(-1:0.001:-0.001,cc2,1)';
meshz(C1,C2,UB')
xlabel('\delta')
ylabel('\delta'' ')
zlabel('Bound')

[minv mini]=min(UB)
n
