%==========================================================================
%
% convergence_analysis  Estimates the order of convergence and the 
% asymptotic error constant for an iterative method.
%
%   [alpha] = convergence_analysis(x_all)
%   [alpha,lambda] = convergence_analysis(x_all)
%   [alpha,lambda,alpha_all,lambda_all] = convergence_analysis(x_all)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-03
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/Root_Finding_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Root_Finding_Methods.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_all       - (n×(n_iter+1) double) all iterates
%
% -------
% OUTPUT:
% -------
%   alpha       - (1×(n_iter+1) double) order of convergence (best
%                 estimate)
%   lambda      - (1×(n_iter+1) double) asymptotic error constant (best
%                 estimate)
%   alpha_all   - (1×(n_iter+1) double) order of convergence estimates at
%                 all iterations
%   lambda_all  - (1×(n_iter+1) double) asymptotic error constant estimates
%                 at all iterations
%
%==========================================================================
function [alpha,lambda,alpha_all,lambda_all] = convergence_analysis(x_all)
    
    % number of iterations
    n_iter = size(x_all,2)-1;
    
    % handles the edge cases where no convergence estimates can be produced
    if n_iter <= 2
        alpha = NaN;
        lambda = NaN;
        alpha_all = NaN(1,n_iter+1);
        lambda_all = NaN(1,n_iter+1);
        return
    end
    
    % preallocates vectors to store order of convergence and asymptotic
    % error constant estimates
    alpha_all = zeros(1,n_iter+1);
    lambda_all = zeros(1,n_iter+1);
    
    % cannot estimate quantites at first two and last iterations
    alpha_all(1) = NaN; lambda_all(1) = NaN;
    alpha_all(2) = NaN; lambda_all(2) = NaN;
    alpha_all(n_iter+1) = NaN; lambda_all(n_iter+1) = NaN;
    
    % estimates order of convergence at all other iterations
    for k = 3:n_iter
        
        % extracts relevant vectors
        a = x_all(:,k+1);
        b = x_all(:,k);
        c = x_all(:,k-1);
        d = x_all(:,k-2);
        
        % calculates residuals
        r_ab = norm(a-b);
        r_bc = norm(b-c);
        r_cd = norm(c-d);
        
        % estimates the order of convergence
        alpha_all(k) = log(r_ab/r_bc)/log(r_bc/r_cd);
        
        % estimates the asymptotic error constant
        lambda_all(k) = r_ab/r_bc^alpha_all(k);
        
    end
    
    % best estimates
    alpha = alpha_all(n_iter);
    lambda = lambda_all(n_iter);
    
end