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
% Last Update: 2023-03-20
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
%   x0      - (1×1 OR 1×2 double) initial guess (x0 = x₀ ∈ ℝ) or initial 
%             interval (x0 = [a,b] ∈ ℝ¹ˣ²)
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • batol      - (1×1 double) absolute bracket tolerance (defaults to
%                      2ε)
%       • max_iter   - (1×1 double) maximimum number solver of iterations
%                      allowed (defaults to 200)
%       • max_feval  - (1×1 double) maximimum number of function
%                      evaluations allowed (defaults to 200)
%       • rebracket  - (1×1 logical) true if initial bracketing interval 
%                      should be updated to ensure sign change, false 
%                      otherwise (defaults to false)
%       • print      - (1×1 logical) true if solver progress should be
%                      printed, false otherwise (defaults to false)
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
%       • x_all      - (1×(n_iter+1) double) root estimates at all 
%                      iterations
%       • a_all      - (1×(n_iter+1) double) bracketing interval lower 
%                      bounds at all iterations
%       • b_all      - (1×(n_iter+1) double) bracketing interval upper 
%                      bounds at at all iterations
%       • n_int_iter - (1×1 double) number of iterations to find a 
%                      bracketing interval
%       • n_iter     - (1×1 double) number of solver iterations
%       • n_feval    - (1×1 double) number of function evaluations
%
%==========================================================================
function [x,output] = root_itp(f,x0,opts)
    
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
    
    % sets absolute bracket tolerance (defaults to 2ε)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'batol')
        batol = 2*eps;
    else
        batol = opts.batol;
    end
    
    % sets maximum number of iterations allowed (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'max_iter')
        max_iter = 200;
    else
        max_iter = opts.max_iter;
    end
    
    % sets maximum number of function evaluations allowed (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'max_feval')
        max_feval = 200;
    else
        max_feval = opts.max_feval;
    end
    
    % determines if the initial bracketing interval should be updated
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'rebracket')
        rebracket = false;
    else
        rebracket = opts.rebracket;
    end
    
    % determines if solver progress should be printed (defaults to no)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'print')
        print = false;
    else
        print = opts.print;
    end
    
    % starts counting function evaluations and number of iterations for
    % finding an initial bracketing interval
    n_feval = 0;
    n_int_iter = 0;
    
    % obtains bracketing interval
    if length(x0) == 1
        [a,b,nf,ni] = bracket_sign_change(f,x0);
        n_feval = n_feval+nf;
        n_int_iter = n_int_iter+ni;
        rebracket = false;
    else
        a = x0(1);
        b = x0(2);
    end
    
    % ensures that a < b
    if (a > b)
        a_old = a;
        a = b;
        b = a_old;
    end
    
    % updates bracketing interval if requested
    if rebracket
        [a,b,nf,ni] = bracket_sign_change(f,[a,b]);
        n_feval = n_feval+nf;
        n_int_iter = n_int_iter+ni;
    end
    
    % auxiliary parameter
    n_max = n0+ceil(log2((b-a)/batol));
    
    % initial guess
    c = (a+b)/2;
    
    % evaluates function at initial guess
    fc = f(c);
    n_feval = n_feval+1;
    
    % returns initial guess if it is a root of f(x) or if the maximum
    % number of function evaluations has already been met/exceeded
    if (fc == 0) || (n_feval >= max_feval)
        x = c;
        output.x_all = c;
        output.a_all = a;
        output.b_all = b;
        output.f_all = fc;
        output.n_int_iter = n_int_iter;
        output.n_iter = 0;
        output.n_feval = n_feval;
        return
    end
    
    % evaluates function at bounds of initial bracketing interval
    ya = f(a);
    yb = f(b);
    n_feval = n_feval+2;
    
    % flips sign on function if f(a) > f(b)
    if ya > yb
        f = @(x) -f(x);
        ya = -ya;
        yb = -yb;
    end
    
    % preallocates arrays to store root estimates, bracketing intervals,
    % and function evaluations
    x_all = zeros(1,max_iter+1);
    a_all = zeros(1,max_iter+1);
    b_all = zeros(1,max_iter+1);
    
    % stores initial guess, bracketing interval, and function evaluation
    x_all(1) = c;
    a_all(1) = a;
    b_all(1) = b;
    
    % prints header for solver progress
    if print
        print_solver_header(true,false);
    end
    
    % iteration
    for k = 1:max_iter
        
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
        r = (batol/2)*2^(n_max-k+1)-(b-a)/2;
        if abs(xt-c) <= r
            x_itp = xt;
        else
            x_itp = c-sigma*r;
        end
        
        % updates interval
        y_itp = f(x_itp);
        n_feval = n_feval+1;
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
        
        % stores kth root estimate, bracketing interval, and function
        % evaluation
        x_all(k+1) = c;
        a_all(k+1) = a;
        b_all(k+1) = b;
        
        % prints solver progress
        if print
            print_solver_progress(k,n_feval,c,[],a,b)
        end
        
        % solver termination on convergence criteria
        if (abs(b-a) <= batol)
            break;
        end
        
        % solver termination on timeout criteria
        if (n_feval >= max_feval)
            break;
        end
        
    end
    
    % prints blank line after last line of solver progress printouts
    if print, fprintf(''); end
    
    % converged root
    x = c;
    
    % number of iterations
    n_iter = k;
    
    % additional outputs
    output.x_all = x_all(1:(n_iter+1));
    output.a_all = a_all(1:(n_iter+1));
    output.b_all = b_all(1:(n_iter+1));
    output.n_int_iter = n_int_iter;
    output.n_iter = n_iter;
    output.n_feval = n_feval;
    
end