%==========================================================================
%
% secant_method  Secant method for finding the root of a univariate, 
% scalar-valued function.
%
%   x = secant_method(f,x0)
%   x = secant_method(f,x0,opts)
%   [x,k] = secant_method(__)
%   [x,k,x_all] = secant_method(__)
%
% See also fzero, bisection_method, secant_method.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-07-06
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
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
%       • k_max      - (1×1 double) maximimum number of iterations 
%                      (defaults to 200)
%       • return_all - (1×1 logical) returns estimates at all iterations if
%                      set to "true"
%       • TOL        - (1×1 double) tolerance (defaults to 10⁻¹⁰)
%
% -------
% OUTPUT:
% -------
%   x       - (1×1 double) root of f(x)
%   k       - (1×1 double) number of solver iterations
%   x_all   - (1×(k+1) double) root estimates at all iterations
%
%==========================================================================
function [x,k,x_all] = secant_method(f,x0,opts)
    
    % ----------------------------------
    % Sets (or defaults) solver options.
    % ----------------------------------
    
    % sets maximum number of iterations (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'k_max')
        k_max = 200;
    else
        k_max = opts.k_max;
    end
    
    % determines if all intermediate estimates should be returned
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'return_all')
        return_all = false;
    else
        return_all = opts.return_all;
    end
    
    % sets tolerance (defaults to 10⁻¹⁰)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'TOL')
        TOL = 1e-10;
    else
        TOL = opts.TOL;
    end
    
    % --------------
    % Secant method.
    % --------------
    
    % root estimates at first and second iterations
    x_prev = x0;
    x_curr = x0+0.001;
    
    % function evaluation at first iteration
    f_prev = f(x0);
    
    % preallocates array and stores estimate at 1st iteration
    if return_all
        x_all = zeros(1,k_max+1);
        x_all(1) = x_prev;
    end
    
    % secant method
    for k = 2:k_max
        
        % stores results in arrays
        if return_all
            x_all(k) = x_curr;
        end
        
        % function evaluation at current iteration
        f_curr = f(x_curr);
        
        % updates root estimate
        x_next = (x_prev*f_curr-x_curr*f_prev)/(f_curr-f_prev);
        
        % terminates solver if converged
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
    
end