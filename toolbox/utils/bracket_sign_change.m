%==========================================================================
%
% bracket_sign_change  Find an interval in which a sign change occurs.
%
%   [a,b] = bracket_sign_change(f,x0)
%   [a,b] = bracket_sign_change(f,[a,b])
%   [a,b] = bracket_sign_change(__,maxiter)
%   [a,b,n_iter,n_feval] = bracket_sign_change(__)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-15
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
%   f           - (1×1 function_handle) univariate, scalar-valued 
%                 function, f(x) (f : ℝ → ℝ)
%   x0          - (1×1 OR 1×2 double) initial guess (x0 = x₀ ∈ ℝ) or 
%                 initial interval (x0 = [a,b] ∈ ℝ¹ˣ²)
%   max_iter    - (OPTIONAL) (1×1 double) maximum number of iterations
%                 allowed (defaults to 200)
%
% -------
% OUTPUT:
% -------
%   a           - (1×1 double) lower bound of interval containing sign
%                 change
%   b           - (1×1 double) lower bound of interval containing sign
%                 change
%   n_iter      - (1×1 double) number of iterations to find a bracketing
%                 interval
%   n_feval     - (1×1 double) number of function evaluations to find a
%                 bracketing interval
%
%==========================================================================
function [a,b,n_iter,n_feval] = bracket_sign_change(f,x0,max_iter)
    
    % sets the initial interval
    if length(x0) == 1
        a = x0;
        b = perturb_iterate(x0);
    else
        a = x0(1);
        b = x0(2);
    end
    
    % returns initial interval if it brackets a sign change
    if f(a)*f(b) < 0
        n_iter = 0;
        n_feval = 2;
        return;
    end
    
    % ensures that a < b
    if (a > b)
        a_old = a;
        a = b;
        b = a_old;
    end
    
    % defaults maximum number of iterations allowed to 200 if not input
    if (nargin < 5) || isempty(max_iter)
        max_iter = 200;
    end
    
    % center of the initial interval
    c = (a+b)/2;
    
    % half-width of the initial interval
    wh = (b-a)/2;
    
    % keeps track of if sign change is found
    sign_change = false;
    
    % keeps expanding interval intil it brackets a sign change
    for k = 1:max_iter
        wh = 2*wh;
        a = c-wh;
        b = c+wh;
        if f(a)*f(b) < 0
            sign_change = true;
            break;
        end
    end
    
    % number of iterations
    n_iter = k;
    
    % number of function evaluations
    n_feval = 2*k+2;
    
    % raises warning if interval with sign change not found
    if ~sign_change
        warning('No interval was found with a sign change.');
    end
    
end