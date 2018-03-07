clc;
clear;
close all;

r=1;
theta=[1];
theta=0.6*theta/(sum(theta));

%Generate Data
nn=10.^[2:0.5:10];

count=0;
a=1;

for n=floor(nn)
 count=count+1;

 %-------------------------------------
 %----Compute exact value -------------
 %-------------------------------------
 noMonteCarlo = 10^6;
 %
 Lam=zeros(1,noMonteCarlo);
 tic
 parfor k=1:noMonteCarlo
%   rng(k); 
  ran=rand(1,n);
  
  %Empirical dist
  q=zeros(1,r+1);
  q(1)=length(find(ran<theta(1)));
  q(2)=n-q(1);
%   for j=2:length(theta)
%    q(j)=length(find(ran<sum(theta(1:j)) & ran>=sum(theta(1:j-1))));
%   end
%   q(r+1)=length(find(ran>=sum(theta(1:r))));
  
  Emp_pmf=q/sum(q);
  
  %Log likelihood ratio
  Lam(k)=2*sum(q.*log((Emp_pmf./[theta,1-sum(theta)])));
 end
 toc
 Ecdf=length(find(Lam<a))/noMonteCarlo;
 Chi_cdf=chi2cdf(a,r);
 
 diff(count)=abs(Ecdf-Chi_cdf);
 
 if(mod(count,1)==0)     
    save(strcat('empirical_5_Bound_',num2str(count),'.mat'));
 end
end

