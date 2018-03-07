clc;
clear;
close all;

r=1;
theta=[1];
theta=0.6*theta/(sum(theta));

%Generate Data
nn=10.^[6.5:0.01:10];
% cn=0;
% for n=floor(nn)
%  cn=cn+1;   
%  Find the upper bound
%  Param
%  
%  c=0;
%  delp=0.1;
%  del=1;
%  tr=1-sum(theta);%Theta_res
%  tm=min(theta,tr); %Theta_min 
%  q=max((r+1)/tm, r*(tm^(-2)-tm^(-1))-1);
%  w=max(2*r, tm^(-2)-tm^(-1) +r*(4+r*( tm^(-2)-tm^(-1) )))/tm;
%  h=400*r^0.25*(1/tm -1)^1.5; 
%  h=400*(r^0.25)*(sum(theta.*(1./theta -1).^1.5) + tr*(1/tr -1)^1.5); 
%  eta=r^5/(tm^6);
%  zeta=r^3/(3*tm^3);
%  dels=delp/(1-delp);
%  a=2;%CHOOSE the worse a later
%  Delp=a*dels + (eta*del^2)/(n*(1-delp))+(zeta*del^1.5)/sqrt(n);
%  
%  eps=2*r*exp(-0.5*n*(delp^2)/(w+delp*q/3));
%  eps2=2*r*exp(-2*del/r);
%  eps3=h/sqrt(n);
%  eps4=Delp*(a/2 + Delp)^(0.5*r-1);
%  UB(cn)=eps+eps2+eps3+eps4;
% end
% loglog(floor(nn),UB)
% hold on;


count=0;
% nn=nn(1:10);
for n=floor(nn)
 count=count+1;

 %-------------------------------------
 %----Compute exact value -------------
 %-------------------------------------
 noMonteCarlo = 10^3;
 %
 Lam=zeros(1,noMonteCarlo);
 tic
 parfor k=1:noMonteCarlo
  rng(k); 
  ran=rand(1,n);
  
  %Empirical dist
  q=zeros(1,r+1);
  q(1)=length(find(ran<theta(1)));
%   for j=2:length(theta)
%    q(j)=length(find(ran<sum(theta(1:j)) & ran>=sum(theta(1:j-1))));
%   end
%   q(r+1)=length(find(ran>=sum(theta(1:r))));
  q(2)=length(find(ran>=sum(theta(1))));
  Emp_pmf=q/sum(q);
  
  %Log likelihood ratio
  LLR_true=log([theta,1-sum(theta)])*q';
  LLR_emp=(log(Emp_pmf))*q';
  Lam(k)=2*(LLR_emp-LLR_true);
 end
 toc
 [Emp_dist x_emp]=ecdf(Lam);
 Lam_chi=chi2cdf(sort(Lam),r);
 
 %Find maximum difference
 x1 = x_emp';
 y1 = Emp_dist';
 x2 = sort(Lam);
 y2 = Lam_chi;
 [x1, index1] = unique(x1);
 [x2, index2] = unique(x2);
 pr_y2 = interp1(x1,y1(index1),x2);
 [maxdist(count),maxidx] = max(abs(pr_y2-y2(index2)));    
 maxdist(count);
 Data{count}=[x2;y2(index2);pr_y2];
 clear('Lam','Rmat');
 if(mod(count,5)==0)     
    save(strcat('empirical_2_Bound_',num2str(count),'.mat'));
 end
end
save('empirical_2_Bound_1-10.mat');
plot(floor(nn),maxdist)
hold on;
% plot(floor(nn),floor(nn).^(-0.5))
title('upper bound wrt n')
xlabel('n')
ylabel('upper bound')
print('Empirical','-r300','-djpeg')
% %PDF plot
% X=linspace(-0.2,20,200);
% [h, f]=hist(Lam,X);
% figure;
% plot(X,h/sum(h))
% hold on;
% % plot(f,h/noMonteCarlo./f)
% Lam_chi1=chi2pdf(sort(Lam),r);
% plot(sort(Lam),Lam_chi1,'k','linewidth',2)
% 
% 



