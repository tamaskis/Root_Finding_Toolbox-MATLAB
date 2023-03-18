%% TEST_convergence.m
% Root-Finding Toolbox
%
% Unit testing of the convergence_analysis function.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-02
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

test_suite = TestSuite('convergence tests',false);



%% EDGE CASES

% -------------
% Edge case #1.
% -------------

% sequence of length 1
x_all = 1;

% performance convergence analysis
[alpha,lambda,alpha_all,lambda_all] = convergence_analysis(x_all);

% unit tests
test_suite.add_test(TestEqual(alpha,nan,'alpha for sequence of length 1'));
test_suite.add_test(TestEqual(lambda,nan,'lambda for sequence of length 1'));
test_suite.add_test(TestEqual(alpha_all,nan,'alpha_all for sequence of length 1'));
test_suite.add_test(TestEqual(lambda_all,nan,'lambda_all for sequence of length 1'));

% -------------
% Edge case #2.
% -------------

% sequence of length 2
x_all = [1,2];

% performance convergence analysis
[alpha,lambda,alpha_all,lambda_all] = convergence_analysis(x_all);

% unit tests
test_suite.add_test(TestEqual(alpha,nan,'alpha for sequence of length 2'));
test_suite.add_test(TestEqual(lambda,nan,'lambda for sequence of length 2'));
test_suite.add_test(TestEqual(alpha_all,[nan,nan],'alpha_all for sequence of length 2'));
test_suite.add_test(TestEqual(lambda_all,[nan,nan],'lambda_all for sequence of length 2'));

% -------------
% Edge case #3.
% -------------

% sequence of length 3
x_all = [1,2,3];

% performance convergence analysis
[alpha,lambda,alpha_all,lambda_all] = convergence_analysis(x_all);

% unit tests
test_suite.add_test(TestEqual(alpha,nan,'alpha for sequence of length 3'));
test_suite.add_test(TestEqual(lambda,nan,'lambda for sequence of length 3'));
test_suite.add_test(TestEqual(alpha_all,[nan,nan,nan],'alpha_all for sequence of length 3'));
test_suite.add_test(TestEqual(lambda_all,[nan,nan,nan],'lambda_all for sequence of length 3'));



%% SUBLINEAR CONVERGENCE

% ------------------
% Sublinear test #1.
% ------------------

% defines the sublinear sequence 1/k for k ∈ [95,100]
k = 95:100;
x_all = 1./k;

% performs convergence analysis
[alpha,lambda,alpha_all,lambda_all] = convergence_analysis(x_all);

% unit tests
test_suite.add_test(TestEqual(alpha,0.9899,'alpha for sublinear sequence #1',4));
test_suite.add_test(TestEqual(lambda,0.8932,'lambda for sublinear sequence #1',4));
test_suite.add_test(TestEqual(alpha_all,[NaN,NaN,0.9897,0.9898,0.9899,NaN],'alpha_all for sublinear sequence #1',4));
test_suite.add_test(TestEqual(lambda_all,[NaN,NaN,0.8915,0.8924,0.8932,NaN],'lambda_all for sublinear sequence #1',4));

% ------------------
% Sublinear test #2.
% ------------------

% defines the sublinear sequence 1/ln(ln(ln(k))) for k ∈ [95,100]
k = 95:100;
x_all = 1./log(log(log(k)));

% performs convergence analysis
[alpha,lambda,alpha_all,lambda_all] = convergence_analysis(x_all);

% unit tests
test_suite.add_test(TestEqual(alpha,0.9871,'alpha for sublinear sequence #2',4));
test_suite.add_test(TestEqual(lambda,0.9206,'lambda for sublinear sequence #2',4));
test_suite.add_test(TestEqual(alpha_all,[NaN,NaN,0.9868,0.9869,0.9871,NaN],'alpha_all for sublinear sequence #2',4));
test_suite.add_test(TestEqual(lambda_all,[NaN,NaN,0.9192,0.9199,0.9206,NaN],'lambda_all for sublinear sequence #2',4));



%% LINEAR CONVERGENCE

% defines the linear sequence 1/2ᵏ for k ∈ [95,100]
k = 95:100;
x_all = 1./2.^k;

% performs convergence analysis
[alpha,lambda,alpha_all,lambda_all] = convergence_analysis(x_all);

% unit tests
test_suite.add_test(TestEqual(alpha,1,'alpha for linear sequence',4));
test_suite.add_test(TestEqual(lambda,0.5,'lambda for linear sequence',4));
test_suite.add_test(TestEqual(alpha_all,[NaN,NaN,1,1,1,NaN],'alpha_all for linear sequence',4));
test_suite.add_test(TestEqual(lambda_all,[NaN,NaN,0.5,0.5,0.5,NaN],'lambda_all for linear sequence',4));



%% QUADRATIC CONVERGENCE

% defines the quadratic sequence xₖ₊₁ = (xₖ)² for k ∈ [25,32] with 
% x₀ = 0.9999999
n_iter = 32;
x_all = zeros(1,n_iter+1);
x_all(1) = 0.9999999;
for k = 1:n_iter
    x_all(k+1)=x_all(k)^2;
end
x_all = x_all(26:33);

% performs convergence analysis
[alpha,lambda,alpha_all,lambda_all] = convergence_analysis(x_all);

% unit tests
test_suite.add_test(TestEqual(alpha,2,'alpha for quadratic sequence',4));
test_suite.add_test(TestEqual(lambda,1,'lambda for quadratic sequence',4));
test_suite.add_test(TestEqual(alpha_all,[NaN,NaN,2.0203,2.0004,2.0000,2.0000,2.0000,NaN],'alpha_all for quadratic sequence',4));
test_suite.add_test(TestEqual(lambda_all,[NaN,NaN,1.1487,1.0049,1.0000,1.0000,1.0000,NaN],'lambda_all for quadratic sequence',4));



%% RUNS TEST SUITE

test_suite.run;