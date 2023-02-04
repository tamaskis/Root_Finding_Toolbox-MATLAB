%==========================================================================
%
% bracket_sign_change  Finds an interval in which a sign change occurs.
%
%   [a,b] = bracket_sign_change(f,x0)
%   [a,b] = bracket_sign_change(f,[a0,b0])
%   [a,b] = bracket_sign_change(__,kappa,maxiter)
%   [a,b,niter,nfeval] = bracket_sign_change(__)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2023-02-04
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
%   x0      - (1×1 OR 1×2 double) initial guess (x0 = x₀ ∈ ℝ) or initial 
%             interval ([a₀,b₀] ∈ ℝ²)
%   kappa   - (OPTIONAL) (1×1 double) scaling factor (defaults to 2)
%   maxiter - (OPTIONAL) (1×1 double) maximum number of iterations
%             (defaults to 200)
%
% -------
% OUTPUT:
% -------
%   a       - (1×1 double) lower bound of interval containing sign change
%   b       - (1×1 double) lower bound of interval containing sign change
%   niter   - (1×1 double) number of iterations to find a bracketing
%             interval
%   nfeval  - (1×1 double) number of function evaluations
%
%==========================================================================
function [a,b,niter,nfeval] = bracket_sign_change(f,x0,kappa,maxiter)
    
    % sets the initial interval
    if length(x0) == 1
        a = x0;
        b = perturb_iterate(x0);
    else
        a = x(1);
        b = x(2);
    end
    
    % returns initial interval if it brackets a sign change
    %   --> no iterations performed
    %   --> 2 function evaluations
    if f(a)*f(b) < 0
        niter = 0;
        nfeval = 2;
        return;
    end
    
    % ensures that a < b
    if (a > b)
        a_old = a;
        a = b;
        b = a_old;
    end
    
    % defaults expansion factor to 2 if not input
    if (nargin < 4) || isempty(kappa)
        kappa = 2;
    end
    
    % defaults maximum number of iterations to 200 if not input
    if (nargin < 5) || isempty(maxiter)
        maxiter = 200;
    end
    
    % center of the initial interval
    c = (a+b)/2;
    
    % half-width of the initial interval
    wh = (b-a)/2;
    
    % keeps expanding interval intil it brackets a sign change
    for k = 1:maxiter
        wh = kappa*wh;
        a = c-wh;
        b = c+wh;
        if f(a)*f(b) < 0
            break;
        end
    end
    
    % number of function evaluations
    nfeval = 2*k+2;
    
    % number of iterations
    niter = k;
    
    % raises warning if interval with sign change not found
    warning('No interval was found with a sign change.');
    
end