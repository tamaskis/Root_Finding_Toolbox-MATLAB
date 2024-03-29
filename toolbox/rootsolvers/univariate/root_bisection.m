%==========================================================================
%
% root_bisection  Bisection method for finding the root of a univariate, 
% scalar-valued function.
%
%   x = root_bisection(f,[a,b])
%   x = root_bisection(f,x0)
%   x = root_bisection(__,opts)
%   [x,output] = root_bisection(__)
%
% See also root_brent_dekker, root_iteration, root_itp, root_newton, 
% root_secant.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-16
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
%       • vtol       - (1×1 double) value tolerance (defaults to 0)
%       • max_iter   - (1×1 double) maximimum number solver of iterations
%                      allowed (defaults to k₁⸝₂ = ⌈log₂((b-a)/batol)⌉)
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
function [x,output] = root_bisection(f,x0,opts)
    
    % sets absolute bracket tolerance (defaults to 2ε)
    %   TODO: do we need this restriction?
    %   • note that we do not allow bracket tolerances below 2ε
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'batol')
        batol = 2*eps;
    else
        batol = max(opts.batol,2*eps);
    end
    
    % sets value tolerance (defaults to 0)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'vtol')
        vtol = 0;
    else
        vtol = opts.vtol;
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
    
    % determines k₁⸝₂
    k_12 = ceil(log2(abs(b-a)/batol));
    
    % sets maximum number of iterations allowed
    if (nargin == 3) && ~isempty(opts) && isfield(opts,'max_iter')
        max_iter = min(opts.max_iter,k_12);
    else
        max_iter = k_12;
    end
    
    % initial guess
    c = (a+b)/2;
    
    % evaluates function at initial guess
    fc = f(c);
    n_feval = n_feval+1;
    
    % returns initial guess if it is a root of f(x) or if the maximum
    % number of function evaluations has already been met/exceeded
    if (abs(fc) <= vtol) || (n_feval >= max_feval)
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
    
    % evaluates function at lower bound of initial bracketing interval
    fa = f(a);
    n_feval = n_feval+1;
    
    % preallocates arrays to store root estimates, bracketing intervals,
    % and function evaluations
    x_all = zeros(1,max_iter+1);
    a_all = zeros(1,max_iter+1);
    b_all = zeros(1,max_iter+1);
    f_all = zeros(1,max_iter+1);
    
    % stores initial guess, bracketing interval, and function evaluation
    x_all(1) = c;
    a_all(1) = a;
    b_all(1) = b;
    f_all(1) = fc;
    
    % prints header for solver progress
    if print
        print_solver_header(true);
    end
    
    % iteration
    for k = 1:max_iter
        
        % updates interval
        if (fa*fc > 0)
            a = c;
            fa = fc;
        else
            b = c;
        end
        
        % updates root estimate
        c = (a+b)/2;
        
        % evaluates function at updated root estimate
        fc = f(c);
        n_feval = n_feval+1;
        
        % stores kth root estimate, bracketing interval, and function
        % evaluation
        x_all(k+1) = c;
        a_all(k+1) = a;
        b_all(k+1) = b;
        f_all(k+1) = fc;
        
        % prints solver progress
        if print
            print_solver_progress(k,n_feval,c,fc,a,b)
        end
        
        % solver termination on convergence criteria
        if (abs(fc) <= vtol)
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
    output.f_all = f_all(1:(n_iter+1));
    output.n_int_iter = n_int_iter;
    output.n_iter = n_iter;
    output.n_feval = n_feval;
    
end