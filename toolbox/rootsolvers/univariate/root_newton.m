%==========================================================================
%
% root_newton  Newton's method for finding the root of a differentiable,
% univariate, scalar-valued function.
%
%   x = root_newton(f,df,x0)
%   x = root_newton(f,df,x0,opts)
%   [x,output] = root_newton(__)
%
% See also root_bisection, root_brent_dekker, root_iteration, root_itp,
% root_secant.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-01-04
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
%   df      - (1×1 function_handle) derivative of f(x) (f' : ℝ → ℝ)
%   x0      - (1×1 double) initial guess for root
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • TOL        - (1×1 double) tolerance (defaults to 10⁻¹⁰)
%       • k_max      - (1×1 double) maximimum number of iterations, kₘₐₓ
%                      (defaults to 200)
%       • return_all - (1×1 logical) returns estimates at all iterations if
%                      set to "true" (defaults to false)
%
% -------
% OUTPUT:
% -------
%   x       - (1×1 double) root of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all   - (1×(k+1) double) root estimates at all iterations
%       • k       - (1×1 double) number of solver iterations
%       • f_count - (1×1 double) number of function evaluations
%       • d_count - (1×1 double) number of derivative evaluations
%
%==========================================================================
function [x,output] = root_newton(f,df,x0,opts)
    
    % sets tolerance (defaults to 10⁻¹⁰)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'TOL')
        TOL = 1e-10;
    else
        TOL = opts.TOL;
    end
    
    % sets maximum number of iterations (defaults to 200)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'k_max')
        k_max = 200;
    else
        k_max = opts.k_max;
    end
    
    % determines if all intermediate estimates should be returned
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'return_all')
        return_all = false;
    else
        return_all = opts.return_all;
    end
    
    % returns initial guess if it is a root of f(x)
    if f(x0) == 0
        x = x0;
        return
    end
    
    % root estimate at first iteration
    x_curr = x0;
    
    % initializes root estimate at next iteration
    x_next = 0;
    
    % preallocates array to store all intermediate solutions
    if return_all
        x_all = zeros(1,k_max+1);
    end
    
    % iteration
    for k = 1:k_max
        
        % stores results
        if return_all
            x_all(k) = x_curr;
        end
        
        % evaluates derivative at current root estimate
        df_curr = df(x_curr);
        
        % perturbs current root estimate if derivative is 0
        if df_curr == 0
            if x_curr ~= 0
                x_curr = x_curr*(1+100*TOL*abs(x_curr));
            else
                x_curr = 100*TOL;
            end
        end
        
        % updates root estimate
        x_next = x_curr-f(x_curr)/df_curr;
        
        % terminates if converged
        if (abs(x_next-x_curr) < TOL)
            break;
        end
        
        % stores updated root estimate for next iteration
        x_curr = x_next;
        
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
    output.f_count = k+1;
    output.d_count = k;
    
end