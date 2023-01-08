%==========================================================================
%
% root_secant  Secant method for finding the root of a univariate, 
% scalar-valued function.
%
%   x = root_secant(f,x0)
%   x = root_secant(f,x0,opts)
%   [x,output] = root_secant(__)
%
% See also root_bisection, root_brent_dekker, root_iteration, root_itp,
% root_newton.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-01-07
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
%   f       - (1×1 function_handle) univariate, scalar-valued function, 
%             f(x) (f : ℝ → ℝ)
%   x0      - (1×1 OR 1×2 double) two options:
%               --> if x₀ ∈ ℝ, then initial guess is input
%               --> if x₀ ∈ ℝ², then initial guess (1st element of x0) and
%                   root estimate at first iteration (2nd element of x0) 
%                   are input
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • TOL        - (1×1 double) tolerance (defaults to 10⁻¹⁰)
%       • k_max      - (1×1 double) maximimum number of iterations, kₘₐₓ
%                      (defaults to 200)
%
% -------
% OUTPUT:
% -------
%   x       - (1×1 double) root of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all   - (1×(k+1) double) root estimates at all iterations
%       • k       - (1×1 double) number of solver iterations
%       • f_count - (1×1 double) number of function evaluations
%
%==========================================================================
function [x,output] = root_secant(f,x0,opts)
    
    % sets tolerance (defaults to 10⁻¹⁰)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'TOL')
        TOL = 1e-10;
    else
        TOL = opts.TOL;
    end
    
    % sets maximum number of iterations (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'k_max')
        k_max = 200;
    else
        k_max = opts.k_max;
    end
    
    % initializes previous and current root estimates
    if length(x0) == 2
        x_prev = x0(1);
        x_curr = x0(2);
    else
        x_prev = x0;
        if x0 ~= 0
            x_curr = x0*(1+100*TOL*abs(x0));
        else
            x_curr = 100*TOL;
        end
    end
    
    % evaluates function at initial guess
    f_prev = f(x_prev);
    
    % returns initial guess if it is a root of f(x)
    if f_prev == 0
        x = x_prev;
        output.x_all = [x_prev,x_curr];
        output.k = 1;
        output.f_count = 1;
        return
    end
    
    % preallocates array to store all intermediate solutions
    x_all = zeros(1,k_max+1);
    
    % stores initial guess and root estimate at 1st iteration
    x_all(1) = x_prev;
    x_all(2) = x_curr;
    
    % counter for number of times current root estimate is perturbed
    n_perturb = 0;
    
    % iteration
    for k = 2:k_max
        
        % evaluates function at current iteration
        f_curr = f(x_curr);
        
        % perturbs current root estimate if function evaluation is same
        if f_curr == f_prev
            if x_curr ~= 0
                x_curr = x_curr*(1+100*TOL*abs(x_curr));
            else
                x_curr = 100*TOL;
            end
            f_curr = f(x_curr);
            n_perturb = n_perturb+1;
        end
        
        % updates root estimate
        x_next = (x_prev*f_curr-x_curr*f_prev)/(f_curr-f_prev);
        
        % stores updated root estimate
        x_all(k+1) = x_next;
        
        % terminates if converged
        if (abs(x_next-x_curr) < TOL)
            break;
        end
        
        % stores next and current root estimates for next iteration
        x_prev = x_curr;
        x_curr = x_next;
        
        % stores current function evaluation for next iteration
        f_prev = f_curr;
        
    end
    
    % converged root
    x = x_next;
    
    % additional outputs
    output.x_all = x_all(1:(k+1));
    output.k = k;
    output.f_count = k+n_perturb+1;
    
end