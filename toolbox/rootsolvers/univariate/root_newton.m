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
%   df      - (1×1 function_handle) derivative of f(x) (f' : ℝ → ℝ)
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
    
    % returns initial guess if it is a root of f(x)
    if f(x0) == 0
        x = x0;
        output.x_all = x;
        output.k = 0;
        output.f_count = 1;
        output.d_count = 0;
        return
    end
    
    % inititalizes current and next root estimates
    x_curr = x0;
    x_next = 0;
    
    % preallocates array to store all intermediate solutions and stores
    % initial guess
    x_all = zeros(1,k_max+1);
    x_all(1) = x0;
    
    % iteration
    for k = 1:k_max
        
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
        
        % stores updated root estimate
        x_all(k+1) = x_next;
        
        % terminates if converged
        if (abs(x_next-x_curr) < TOL)
            break;
        end
        
        % stores updated root estimate for next iteration
        x_curr = x_next;
        
    end
    
    % converged root
    x = x_next;
    
    % additional outputs
    output.x_all = x_all(1:(k+1));
    output.k = k;
    output.f_count = k+1;
    output.d_count = k;
    
end