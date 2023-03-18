%% TEST_root_bisection.m
% Root-Finding Toolbox
%
% Unit testing of the root_bisection function.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-18
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% DEPENDENCIES:
%   • Simple Unit Testing Toolbox (https://tamaskis.github.io/Simple_Unit_Testing_Toolbox-MATLAB)
%   • Root-Finding Toolbox (https://tamaskis.github.io/Root_Finding_Toolbox-MATLAB/)



%% SCRIPT SETUP

% clears Workspace and closes all figures
clear; close all;



%% INITIALIZE TEST SUITE

test_suite = TestSuite('root_bisection tests',false);



%% ROOT AT MIDPOINT OF INITIAL INTERVAL

% function with root at x = 1
f = @(x) x^2-1;

% initial interval with midpoint at 1
x0 = [-9,11];

% solve for root
[x,output] = root_bisection(f,x0);

% unit tests
test_suite.add_test(TestEqual(x,1,'initial interval midpoint = root --> root'));
test_suite.add_test(TestEqual(output.n_iter,0,'initial interval midpoint = root --> number of iterations'));
test_suite.add_test(TestEqual(output.n_feval,1,'initial interval midpoint = root --> number of function evaluations'));



%% COUNTING FUNCTION EVALUTIONS

% function with root at x = 1
f = @(x) x^2-1;

% initial interval
x0 = [0,9999999];

% function handle for function that returns n_feval
g = @(f) return_n_feval(f,x0);

% unit test
test_suite.add_test(TestFunctionCount(f,g,...
    'counting function evaluations'));



%% TESTING TOLERANCE


%% TEST MAX ITERATIONS


%% TEST MAX F EVAL



%% RUNS TEST SUITE

test_suite.run;



%% AUXILIARY FUNCTIONS

function n_feval = return_n_feval(f,x0)
    [~,output] = root_bisection(f,x0);
    n_feval = output.n_feval;
end