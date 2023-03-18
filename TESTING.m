% TODO:
%   --> https://www.mathworks.com/help/matlab/ref/try.html
%   --> https://www.mathworks.com/matlabcentral/answers/142960-using-try-catch-to-get-warning-message
%   --> https://www.mathworks.com/matlabcentral/answers/364719-detect-warning-and-take-action
%   --> https://www.mathworks.com/matlabcentral/answers/419823-how-to-catch-warnings

clear;clc;close all;


%%
% s = warning('error','MATLAB:singularMatrix');
% A = zeros(4);
% b = [1;2;3;4];
% 
% 
% try
%     x = 1;
%     A\b
% catch ME
%     x
%     ME.identifier
% end
% 
% warning(s)
% 
% A\b








%%

addpath(genpath('toolbox'));

f = @(x) x^2-1;
df = @(x) 2*x;
opts.print = false;
[x,output] = root_bisection(f,[0,1000],opts);
%[x,output] = root_newton(f,df,10,opts);

[alpha_all,lambda_all] = convergence_analysis(output.x_all)

convergence_plot(output.x_all)