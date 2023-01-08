clear; clc; close all;

% f = @(x) exp(-exp(-x))-x;
% x0 = [0,1];
% 
% f_matlab = @(f,x0) fzero(f,x0);
% 
% f_rft = @(f,x0) root_brent_dekker(f,x0);
% 
% TIME_EVALUATION(f_matlab,{f,x0});
% TIME_EVALUATION(f_rft,{f,x0});
% 
% 
% x1 = f_matlab(f,x0)
% 
% x2 = f_rft(f,x0)
% 
% x1-x2


fun = @(x) exp(-exp(-x)) - x; % function
x0 = [0 1]; % initial interval
options = optimset('Display','iter'); % show iterations
[x fval exitflag output] = fzero(fun,x0,options)

opts.return_all = true;
[x,output] = root_brent_dekker(fun,x0,opts)