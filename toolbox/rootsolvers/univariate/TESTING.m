clear; clc; close all;

f = @(x) x^2-1;
x0 = 1000;

df = @(x) 2*x;

opts.return_all = true;
f_fast = @(f,df,x0) root_newton(f,df,x0);
f_slow = @(f,df,x0) root_newton(f,df,x0,opts);

TIME_EVALUATION(f_fast,{f,df,x0});
TIME_EVALUATION(f_slow,{f,df,x0});