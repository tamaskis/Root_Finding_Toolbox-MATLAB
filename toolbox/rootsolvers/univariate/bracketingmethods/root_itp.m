%==========================================================================
%
% root_itp  Interpolate, Truncate, and Project method for finding the root 
% of a univariate, scalar-valued function.
%
%   x = root_itp(f,[a,b])
%   x = root_itp(f,x0)
%   x = root_itp(__,opts)
%   [x,output] = root_itp(__)
%
% See also root_bisection, root_brent_dekker, root_iteration, root_newton, 
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
%       • rebracket  - (1×1 double) true if initial bracket should be
%                      updated to ensure sign change, false otherwise 
%                      (defaults to false)
%       • kappa1     - (1×1 double) tuning parameter, κ₁ (defaults to 0.1)
%                       --> NOTE: κ₁ ∈ (0,∞)
%       • kappa2     - (1×1 double) tuning parameter, κ₂ (defaults to 
%                      0.98(1+φ), where φ = (1+√(5))/2 is the golden ratio)
%                       --> NOTE: κ₂ ∈ [1,1+φ)
%       • n0         - (1×1 double) tuning parameter, n₀ (defaults to 1)
%                       --> NOTE: n₀ ∈ [0,∞)
%
% -------
% OUTPUT:
% -------
%   x       - (1×1 double) root of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all   - (1×(k+1) double) root estimates at all iterations
%       • a_all   - (1×(k+1) double) bracketing interval lower bounds at 
%                   all iterations
%       • b_all   - (1×(k+1) double) bracketing interval upper bounds at 
%                   all iterations
%       • k       - (1×1 double) number of solver iterations
%       • f_count - (1×1 double) number of function evaluations
%
%==========================================================================
function [x,output] = root_itp(f,x0,opts)
    
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
    
    % determines if the initial bracketing interval should be updated
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'rebracket')
        rebracket = false;
    else
        rebracket = opts.rebracket;
    end
    
    % sets κ₁ (defaults to 0.1)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'kappa1')
        kappa1 = 0.1;
    else
        kappa1 = opts.kappa1;
    end
    
    % sets κ₂ (defaults to 0.98[1+(1+√(5))/2])
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'kappa2')
        kappa2 = 0.98*(1+(1+sqrt(5))/2);
    else
        kappa2 = opts.kappa2;
    end
    
    % sets n₀ (defaults to 1)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'n0')
        n0 = 1;
    else
        n0 = opts.n0;
    end
    
    % obtains bracketing interval
    if length(x0) == 1
        [a,b,f_count_1] = bracket_sign_change(f,x0);
        rebracket = false;
    else
        a = x0(1);
        b = x0(2);
        f_count_1 = 0;
    end
    
    % updates bracketing interval if requested
    if rebracket
        [a,b,f_count_2] = bracket_sign_change(f,[a,b]);
    else
        f_count_2 = 0;
    end
    
    % root estimate at first iteration
    c = (a+b)/2;
    
    % returns root estimate at first iteration if it is a root of f(x)
    if f(c) == 0
        x = c;
        output.x_all = x;
        output.a_all = a;
        output.b_all = b;
        output.k = 0;
        output.f_count = f_count_1+f_count_2+1;
        return
    end
    
    % preallocates arrays to store all intermediate solutions and
    % bracketing intervals
    x_all = zeros(1,k_max+1);
    a_all = zeros(1,k_max+1);
    b_all = zeros(1,k_max+1);
    
    % stores initial guess and bracketing interval
    x_all(1) = x0;
    a_all(1) = a;
    b_all(1) = b;
    
    % evaluates function at bounds of initial bracketing interval
    ya = f(a);
    yb = f(b);
    
    % flips sign on function if f(a) > f(b)
    if ya > yb
        f = @(x) -f(x);
        ya = -ya;
        yb = -yb;
    end
    
    % auxiliary parameter
    n_max = n0+ceil(log2((b-a)/TOL));
    
    % iteration
    for k = 1:k_max
        
        % interpolation
        xf = (yb*a-ya*b)/(yb-ya);
        
        % truncation
        sigma = sign(c-xf);
        delta = kappa1*(b-a)^kappa2;
        if delta <= abs(c-xf)
            xt = xf+sigma*delta;
        else
            xt = c;
        end
        
        % projection
        r = (TOL/2)*2^(n_max-k+1)-(b-a)/2;
        if abs(xt-c) <= r
            x_itp = xt;
        else
            x_itp = c-sigma*r;
        end
        
        % updates interval
        y_itp = f(x_itp);
        if y_itp > 0
            b = x_itp;
            yb = y_itp;
        elseif y_itp < 0
            a = x_itp;
            ya = y_itp;
        else
            a = x_itp;
            b = x_itp;
        end
        
        % updates root estimate
        c = (a+b)/2;
        
        % stores updated root estimate and bracketing interval
        x_all(k+1) = c;
        a_all(k+1) = a;
        b_all(k+1) = b;
        
        % terminates if converged
        if ((b-a) < TOL)
            break;
        end
        
    end
    
    % converged root
    x = c;
    
    % additional outputs
    output.x_all = x_all(1:(k+1));
    output.a_all = a_all(1:(k+1));
    output.b_all = b_all(1:(k+1));
    output.k = k;
    output.f_count = f_count_1+f_count_2+3+k;
    
end