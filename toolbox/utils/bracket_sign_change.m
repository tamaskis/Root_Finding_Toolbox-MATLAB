%==========================================================================
%
% bracket_sign_change  Finds an interval in which a sign change occurs.
%
%   [a,b] = bracket_sign_change(f,x0)
%   [a,b] = bracket_sign_change(f,[a0,b0])
%   [a,b] = bracket_sign_change(__,kappa,k_max)
%   [a,b,f_count] = bracket_sign_change(__)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2023-01-05
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
%   f       - (1×1 function_handle) univariate function, f(x) (f : ℝ → ℝ)
%   x0      - (1×1 OR 1×2 double) two options:
%               --> if x₀ ∈ ℝ, then default to [a,b] = [x₀,x₀+1]
%               --> if x₀ ∈ ℝ², then a = x₁ and b = x₂
%   kappa   - (OPTIONAL) (1×1 double) scaling factor (defaults to 2)
%   k_max   - (OPTIONAL) (1×1 double) maximum number of expansions
%             (defaults to 200)
%
% -------
% OUTPUT:
% -------
%   a       - (1×1 double) lower bound of interval containing sign change
%   b       - (1×1 double) lower bound of interval containing sign change
%   f_count - (1×1 double) number of function evaluations
%
%==========================================================================
function [a,b,f_count] = bracket_sign_change(f,x0,kappa,k_max)
    
    % sets the initial interval
    if length(x0) == 2
        a = x(1);
        b = x(2);
    else
        a = x0;
        b = x0+1;
    end
    
    % returns initial interval if it brackets a sign change
    if f(a)*f(b) < 0
        return;
    end
    
    % ensures that a < b
    if (a > b)
        a_old = a;
        a = b;
        b = a_old;
    end
    
    % defaults expansion factor to 2 if not input
    if (nargin < 3) || isempty(kappa)
        kappa = 2;
    end
    
    % defaults maximum number of iterations to 200 if not input
    if (nargin < 4) || isempty(k_max)
        k_max = 200;
    end
    
    % center of the initial interval
    c = (a+b)/2;
    
    % half-width of the initial interval
    wh = (b-a)/2;
    
    % keeps expanding interval intil it brackets a sign change
    for k = 1:k_max
        wh = kappa*wh;
        a = c-wh;
        b = c+wh;
        if f(a)*f(b) < 0
            f_count = 2+2*k;
            return;
        end
    end
    
    % raise warning if interval with sign change not found
    warning('No interval was found with a sign change.');
    
end