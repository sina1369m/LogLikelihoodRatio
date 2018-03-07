clear;
clc;
close all;

load('UpperBound4.mat', 'UB','UB_','c','nn' );
N_ub=floor(nn);
clear nn;

load('empirical_4_Bound_332.mat', 'diff','nn');
N_emp=floor(nn(1:length(diff)));
clear nn;

s=find(UB>1);
last_n=N_ub(max(s));
s_=find(UB_>1);

t=find(N_emp>last_n);



loglog(N_ub(max(s):end),UB(max(s):end),'-r');
hold on;
loglog(N_ub(max(s_):end),UB_(max(s_):end),'--r');
%Trivial bound
loglog([N_emp(1),N_emp(min(t))],[1 1],'-r')
loglog(N_emp,diff,'-b')


grid on
xlabel('n')
ylabel('The bound')
legend('Upper bound Theorem 1',strcat('Refined bound with c=',num2str(c)),'Trivial Bound','Simulated difference')
title('a=1')
print('Comp_Bound','-r300','-djpeg')


