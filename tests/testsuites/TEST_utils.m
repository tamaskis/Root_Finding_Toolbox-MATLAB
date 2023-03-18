%% TEST_utils.m
% Root-Finding Toolbox
%
% Unit testing of the bracket_sign_change and perturb_iterate functions.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-03
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

test_suite = TestSuite('utils tests',false);



%% perturb_iterate

test_suite.add_test(TestEqual(perturb_iterate(0),100*eps,'perturb x = 0'));
test_suite.add_test(TestEqual(perturb_iterate(1),1.000000000000022,'perturb x = 1',15));
test_suite.add_test(TestEqual(perturb_iterate(99999),99999.00000222040,'perturb x = 1',11));
test_suite.add_test(TestEqual(perturb_iterate([0;0]),[100*eps;100*eps],'perturb x = (0,0)ᵀ'));
test_suite.add_test(TestEqual(perturb_iterate([100;100]),[100.0000000000031;100.0000000000031],'perturb x = (0,0)ᵀ',13));
test_suite.add_test(TestEqual(perturb_iterate([1e5;5e3]),[100000.0000022232;5000.0000001112],'perturb x = (0,0)ᵀ',10));



%% bracket_sign_change

% test #1
f = @(x) x;
x0 = [-1,1];
[a,b,n_iter,n_feval] = bracket_sign_change(f,x0);
test_suite.add_test(TestEqual([a,b],[-1,1],'interval for f(x) = x and x₀ = [-1,1]',4));
test_suite.add_test(TestEqual(n_iter,0,'n_iter for f(x) = x and x₀ = [-1,1]'));
test_suite.add_test(TestEqual(n_feval,2,'n_feval for f(x) = x and x₀ = [-1,1]'));

% test #2
f = @(x) x;
x0 = 0;
[a,b,n_iter,n_feval] = bracket_sign_change(f,x0);
test_suite.add_test(TestEqual([a,b],[-0.1110e-13,0.3331e-13],'interval for f(x) = x and x₀ = 0'));
test_suite.add_test(TestEqual(n_iter,1,'n_iter for f(x) = x and x₀ = 0'));
test_suite.add_test(TestEqual(n_feval,4,'n_feval for f(x) = x and x₀ = 0'));

% test #3
f = @(x) x;
x0 = 10;
[a,b,n_iter,n_feval] = bracket_sign_change(f,x0);
test_suite.add_test(TestEqual([a,b],[-5.6250,25.6250],'interval for f(x) = x and x₀ = 10',4));
test_suite.add_test(TestEqual(n_iter,47,'n_iter for f(x) = x and x₀ = 10'));
test_suite.add_test(TestEqual(n_feval,96,'n_feval for f(x) = x and x₀ = 10'));

% test #4
f = @(x) x;
x0 = [50,100];
[a,b,n_iter,n_feval] = bracket_sign_change(f,x0);
test_suite.add_test(TestEqual([a,b],[-25,175],'interval for f(x) = x and x₀ = 10'));
test_suite.add_test(TestEqual(n_iter,2,'n_iter for f(x) = x and x₀ = 10'));
test_suite.add_test(TestEqual(n_feval,6,'n_feval for f(x) = x and x₀ = 10'));



%% RUNS TEST SUITE

test_suite.run;