clear;
clc;

load('Epirical_Bound_1-4.mat');

n_T=nn;
Emp_Bound=maxdist;

clearvars -except n_T Emp_Bound;
%%
load('Epirical_Bound_5.mat');
n_T=[n_T nn];
Emp_Bound=[Emp_Bound maxdist];

clearvars -except n_T Emp_Bound;
%%
plot(floor(n_T),Emp_Bound)
