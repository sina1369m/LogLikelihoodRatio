clc;
clear;
close all;

r=1;
theta=[4];
theta=0.9*theta/(sum(theta));

%Generate Data
n=100;

noMonteCarlo = 5000*n;


    
for k=1:noMonteCarlo
    rng(k);
    ran=rand(1,n);    
    q(1)=length(find(ran<theta(1)));
    for j=2:length(theta)
        q(j)=length(find(ran<sum(theta(1:j)) & ran>=sum(theta(1:j-1))));
    end
    q(r+1)=length(find(ran>=sum(theta(1:r))));
    
    Emp_pmf=q/sum(q);
    LLR_true=log([theta,1-sum(theta)])*q';
    LLR_emp=(log(Emp_pmf))*q';
    Lam(k)=2*(LLR_emp-LLR_true);
end

%CDF plot
figure;
ecdf(Lam);
h1=cdfplot(Lam);
hold on;

Lam_chi=chi2cdf(sort(Lam),r);
h2=plot(sort(Lam),Lam_chi,'k','linewidth',1.2);

%Find maximum difference
x1 = get(h1, 'XData');
y1 = get(h1, 'YData');
x2 = get(h2, 'XData');
y2 = get(h2, 'YData');

x1=x1(2:end-1);
x2=x2(2:end-1);
y1=y1(2:end-1);
y2=y2(2:end-1);

[x1, index1] = unique(x1);
[x2, index2] = unique(x2);

pr_y2 = interp1(x1,y1(index1),x2);
[maxdist,maxidx] = max(abs(pr_y2-y2(index2)));    
maxdist;

xlabel('\lambda')
ylabel('F(\lambda)')
text(12,0.15,strcat('Maximum difference=',num2str(maxdist)))
print('CompCDF-n100','-r300','-djpeg')

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



