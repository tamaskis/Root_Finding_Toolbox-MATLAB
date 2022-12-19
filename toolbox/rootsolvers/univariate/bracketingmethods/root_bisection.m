%==========================================================================
%
% root_bisection  Bisection method for finding the root of a univariate, 
% scalar-valued function.
%
%   x = root_bisection(f,[a,b])
%   x = root_bisection(f,x0)
%   x = root_bisection(__,opts)
%   [x,k] = root_bisection(__)
%   [x,k,x_all] = root_bisection(__)
%
% See also root_brent_dekker, root_iteration, root_itp, root_newton, 
% root_secant.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-12-11
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
%               --> if x₀ ∈ ℝ (i.e. initial guess input), then we attempt
%                   to use bracket_sign_change to find an initial 
%                   bracketing interval
%               --> if x₀ ∈ ℝ² (i.e. initial bracketing interval input), 
%                   then a = x₁ and b = x₂
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • TOL        - (1×1 double) tolerance (defaults to 10⁻¹⁰)
%       • k_max      - (1×1 double) maximimum number of iterations, kₘₐₓ
%                      (defaults to 200)
%       • return_all - (1×1 logical) returns estimates at all iterations if
%                      set to "true" (defaults to false)
%       • rebracket  - (1×1 double) true if initial bracket should be
%                      updated to ensure sign change, false otherwise 
%                      (defaults to false)
%
% -------
% OUTPUT:
% -------
%   x       - (1×1 double) root of f(x)
%   k       - (1×1 double) number of solver iterations
%   x_all   - (1×(k+1) double) root estimates at all iterations
%
%==========================================================================
function [x,k,x_all] = root_bisection(f,x0,opts)
    
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
    
    % determines if all intermediate estimates should be returned
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'return_all')
        return_all = false;
    else
        return_all = opts.return_all;
    end
    
    % determines if the initial bracketing interval should be updated
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'rebracket')
        rebracket = false;
    else
        rebracket = opts.rebracket;
    end
    
    % obtains bracketing interval
    if length(x0) == 1
        [a,b] = bracket_sign_change(f,x0);
        rebracket = false;
    else
        a = x0(1);
        b = x0(2);
    end
    
    % updates bracketing interval if requested
    if rebracket
        [a,b] = bracket_sign_change(f,[a,b]);
    end
    
    % root estimate at first iteration
    c = (a+b)/2;
    
    % returns root estimate at first iteration if it is a root of f(x)
    if f(c) == 0
        x = c;
        return
    end
    
    % function evaluations at first iteration
    fa = f(a);
    fc = f(c);
    
    % preallocates array to store all intermediate solutions
    if return_all
        x_all = zeros(1,k_max+1);
    end
    
    % iteration
    for k = 1:k_max
        
        % stores results
        if return_all
            x_all(k) = c;
        end
        
        % updates interval
        if fc == 0
            break;
        elseif (fa*fc > 0)
            a = c;
            fa = fc;
        else
            b = c;
        end
        
        % updates root estimate
        c = (a+b)/2;
        
        % terminates if converged
        if ((b-a) < TOL)
            break;
        end
        
        % function evaluation at updated root estimate
        fc = f(c);
        
    end
    
    % converged root
    x = c;
    
    % stores converged result and trims array
    if return_all
        x_all(k+1) = x;
        x_all = x_all(1:(k+1));
    end
    
end