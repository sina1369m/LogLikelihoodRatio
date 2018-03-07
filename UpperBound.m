clear;
clc;
close all


r=1;
theta=[1];
theta=0.4*theta/(sum(theta));

%Generate Data
nn=10.^[3:0.02:10];
cn=0;
for n=floor(nn)
 cn=cn+1;   
 %Find the upper bound
 %Param
 ccc=0;
 for c1=-0.5:0.001:-0.001
    for c2=0.001:0.01:0.99
     ccc=ccc+1;   
%      c1=-0.2;
%      c2=0.2;
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

     eps=2*r*exp(-0.5*n*(delp^2)/(w+delp*q/3));
     eps2=2*r*exp(-2*del/r);
     eps3=h/sqrt(n);
     c=10;
     eps3_=10/n;
     eps4=Delp*(a/2 + Delp)^(0.5*r-1);
     b_temp(ccc)=eps+eps2+eps3+eps4;   
     b_temp_(ccc)=eps+eps2+eps3_+eps4;   
    end
 end
 [UB(cn) m_index(cn)]=min(b_temp);
 [UB_(cn) m_index_(cn)]=min(b_temp_);
clear b_temp b_temp_;
end

save('UpperBound4.mat')
loglog(floor(nn),UB)
grid on
xlabel('n')
ylabel('The bound')
print('Bound4','-r300','-djpeg')

figure;
loglog(floor(nn),UB)
hold on;
loglog(floor(nn),UB_)
grid on
xlabel('n')
ylabel('Upper bounds')
legend('Upper bound Theorem 1',strcat('Refined bound with c=',num2str(c)))
print('Bound4_2','-r300','-djpeg')


