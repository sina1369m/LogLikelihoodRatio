clear;
clc;
close all


r=1;
theta=[1];
theta=0.6*theta/(sum(theta));
tr=1-sum(theta);
tm=min(min(theta),tr);


n=10.^[5:0.1:10];

h=400*(r^0.25)*(sum(theta.*(1./theta -1).^1.5) + tr*(1/tr -1)^1.5); 
loglog(n,h*n.^(-0.5));
grid on;

title('S3 wrt n')
xlabel('n')
ylabel('S3')
print('S3','-r300','-djpeg')