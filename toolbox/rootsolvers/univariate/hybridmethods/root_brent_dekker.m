%==========================================================================
%
% root_brent_dekker  Brent-Dekker method for finding the root of a 
% univariate, scalar-valued function.
%
%   x = root_brent_dekker(f,[a,b])
%   x = root_brent_dekker(f,x0)
%   x = root_brent_dekker(__,opts)
%   [x,k] = root_brent_dekker(__)
%   [x,k,x_all] = root_brent_dekker(__)
%
% See also root_brent_dekker, root_iteration, root_itp, root_newton, 
% root_secant.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-12-24
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
function [x,k,x_all] = root_brent_dekker(f,x0,opts)
    
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
    
    % returns midpoint of bracketing interval if it is a root of f(x)
    if f((a+b)/2) == 0
        x = (a+b)/2;
        return
    end
    
    % function evaluations at left and right endpoints of bracketing 
    % interval
    fa = f(a);
    fb = f(b);
    
    % initializes function evaluation at previous root estimate
    fc = fb;
    
    % preallocates array to store all intermediate solutions
    if return_all
        x_all = zeros(1,k_max+1);
    end
    
    % iteration
    for k = 1:k_max
        
        % stores results
        if return_all
            x_all(k) = b;
        end
        
        % auxiliary parameters
        if (k == 1) || ((fb > 0) == (fc > 0))
            c = a;
            fc = fa;
            d = b-a;
            e = d;
        end
        
        % updates interval
        if abs(fc) < abs(fb)
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        end
        
        % algorithm tolerance
        delta = 2*eps*abs(b)+TOL;
        
        % auxiliary parameter used for termination, bisection, and 
        % interpolation
        m = 0.5*(c-b);
        
        % terminates if converged
        if (abs(m) < delta) || (fb == 0)
            break;
        end
        
        % performs bisection
        if (abs(e) < delta) || (abs(fa) < abs(fb))
            d = m;
            e = m;
            
        % performs interpolation
        else
            
            % auxiliary parameter for interpolation
            s = fb/fa;
            
            % linear interpolation
            if a == c
                p = 2*m*s;
                q = 1-s;
                
            % inverse quadratic interpolation
            else
                q = fa/fc;
                r = fb/fc;
                p = s*(2*m*q*(q-r)-(b-a)*(r-1));
                q = (q-1)*(r-1)*(s-1);
                
            end
            
            % adjusts sign of p and q
            if p > 0
                q = -q;
            else
                p = -p;
            end
            
        end

        % determines whether to keep results of interpolation or to revert
        % to performing bisection instead
        s = e;
        e = d;
        if (2*p < (3*m*q-abs(delta*q))) || (p < abs(0.5*s*q))
            d = p/q;
        else
            d = m;
            e = m;
        end
        
        % updates root estimate
        a = b;
        fa = fb;
        if abs(d) > delta
            b = b+d;
        elseif m > 0
            b = b+delta;
        else
            b = b-delta;
        end
        
        % function evaluation at updated root estimate
        fb = f(b);
        
    end
    
    % converged root
    x = b;
    
    % stores converged result and trims array
    if return_all
        x_all(k+1) = x;
        x_all = x_all(1:(k+1));
    end
    
end