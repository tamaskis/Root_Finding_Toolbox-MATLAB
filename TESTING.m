% TODO:
%   --> https://www.mathworks.com/help/matlab/ref/try.html
%   --> https://www.mathworks.com/matlabcentral/answers/142960-using-try-catch-to-get-warning-message
%   --> https://www.mathworks.com/matlabcentral/answers/364719-detect-warning-and-take-action
%   --> https://www.mathworks.com/matlabcentral/answers/419823-how-to-catch-warnings

clear;clc;close all;

s = warning('error','MATLAB:singularMatrix');
A = zeros(4);
b = [1;2;3;4];


try
    x = 1;
    A\b
catch ME
    x
    ME.identifier
end

warning(s)

A\b