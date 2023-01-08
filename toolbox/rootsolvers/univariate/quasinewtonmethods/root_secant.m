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
%   x0      - (1×1 double) initial guess for root
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
    
    % function evaluation at first iteration
    f_prev = f(x0);
    
    % returns initial guess if it is a root of f(x)
    if f_prev == 0
        x = x0;
        return
    end
    
    % root estimates at first and second iterations
    x_prev = x0;
    if x0 ~= 0
        x_curr = x0*(1+100*TOL*abs(x0));
    else
        x_curr = 100*TOL;
    end
    
    % preallocates array to store all intermediate solutions and stores 
    % estimate at 1st iteration
    if return_all
        x_all = zeros(1,k_max+1);
        x_all(1) = x_prev;
    end
    
    % counter for number of times current root estimate is perturbed
    n_perturb = 0;
    
    % iteration
    for k = 2:k_max
        
        % stores results
        if return_all
            x_all(k) = x_curr;
        end
        
        % function evaluation at current iteration
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
    
    % stores converged result and trims array
    if return_all
        x_all(k+1) = x;
        x_all = x_all(1:(k+1));
    end
    
    % output structure
    if return_all, output.x_all = x_all; end
    output.k = k;
    output.f_count = k+n_perturb;
    
end