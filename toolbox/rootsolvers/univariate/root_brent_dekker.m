%==========================================================================
%
% root_brent_dekker  Brent-Dekker method for finding the root of a 
% univariate, scalar-valued function.
%
%   x = root_brent_dekker(f,[a,b])
%   x = root_brent_dekker(f,x0)
%   x = root_brent_dekker(__,opts)
%   [x,output] = root_brent_dekker(__)
%
% See also root_brent_dekker, root_iteration, root_itp, root_newton, 
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
%       • tol        - (1×1 double) tolerance (defaults to ε)
%       • vtol       - (1×1 double) value tolerance (defaults to 0)
%       • max_iter   - (1×1 double) maximimum number solver of iterations
%                      allowed (defaults to 200)
%       • max_feval  - (1×1 double) maximimum number of function
%                      evaluations allowed (defaults to 200)
%       • rebracket  - (1×1 logical) true if initial bracketing interval 
%                      should be updated to ensure sign change, false 
%                      otherwise (defaults to false)
%       • print      - (1×1 logical) true if solver progress should be
%                      printed, false otherwise (defaults to false)
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
%       • f_all      - (1×(n_iter+1) double) function evaluations at all
%                      iterations
%       • n_int_iter - (1×1 double) number of iterations to find a 
%                      bracketing interval
%       • n_iter     - (1×1 double) number of solver iterations
%       • n_feval    - (1×1 double) number of function evaluations
%
%==========================================================================
function [x,output] = root_brent_dekker(f,x0,opts)
    
    % sets tolerance (defaults to ε)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'tol')
        tol = 2*eps;
    else
        tol = opts.tol;
    end
    
    % sets value tolerance (defaults to 0)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'vtol')
        vtol = 0;
    else
        vtol = opts.vtol;
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
    % (defaults to no)
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
    
    % evaluates function at upper bound of initial bracketing interval
    fb = f(b);
    
    % returns upper bound of initial bracketing interval if it is a root of
    % f(x) or if the maximum number of function evaluations has already 
    % been met/exceeded
    if (abs(fb) <= vtol) || (n_feval >= max_feval)
        x = b;
        output.x_all = b;
        output.a_all = a;
        output.b_all = b;
        output.f_all = fb;
        output.n_int_iter = n_int_iter;
        output.n_iter = 0;
        output.n_feval = n_feval;
        return
    end
    
    % evaluates function at lower bound of initial bracketing interval
    fa = f(a);
    n_feval = n_feval+1;
    
    % initializes function evaluation at previous root estimate
    fc = fb;
    
    % preallocates arrays to store root estimates, bracketing intervals,
    % and function evaluations
    x_all = zeros(1,max_iter+1);
    a_all = zeros(1,max_iter+1);
    b_all = zeros(1,max_iter+1);
    f_all = zeros(1,max_iter+1);
    
    % stores initial guess, bracketing interval, and function evaluation
    x_all(1) = b;
    a_all(1) = a;
    b_all(1) = b;
    f_all(1) = fb;
    
    % prints header for solver progress
    if print
        print_solver_header(true);
    end
    
    % iteration
    for k = 1:max_iter
        
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
        delta = 2*eps*abs(b)+tol;
        
        % auxiliary parameter used for termination, bisection, and 
        % interpolation
        m = 0.5*(c-b);
        
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
        
        % evaluates function at updated root estimate
        fb = f(b);
        n_feval = n_feval+1;
        
        % stores kth root estimate, bracketing interval, and function
        % evaluation
        x_all(k+1) = b;
        a_all(k+1) = a;
        b_all(k+1) = b;
        f_all(k+1) = fb;
        
        % prints solver progress
        if print
            print_solver_progress(k,n_feval,b,fb,a,b)
        end
        
        % solver termination on convergence criteria
        if (abs(fb) <= vtol) || (abs(m) <= delta)
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
    x = b;
    
    % number of iterations
    n_iter = k;
    
    % additional outputs
    output.x_all = x_all(1:(n_iter+1));
    output.a_all = a_all(1:(n_iter+1));
    output.b_all = b_all(1:(n_iter+1));
    output.f_all = f_all(1:(n_iter+1));
    output.n_int_iter = n_int_iter;
    output.n_iter = n_iter;
    output.n_feval = n_feval;
    
end