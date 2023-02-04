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
%   f       - (1×1 function_handle) univariate, scalar-valued function, 
%             f(x) (f : ℝ → ℝ)
%   x0      - (1×1 OR 1×2 double) initial guess (x0 = x₀ ∈ ℝ) or initial 
%             interval ([a₀,b₀] ∈ ℝ²)
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • xatol      - (1×1 double) absolute step tolerance
%                      (defaults to 10⁻¹⁰)
%       • maxiter    - (1×1 double) maximimum number of iterations 
%                      (defaults to 200)
%       • rebracket  - (1×1 logical) true if initial bracketing interval 
%                      should be updated to ensure sign change, false 
%                      otherwise (defaults to false)
%       • print      - (1×1 logical) true if solver progress should be
%                      printed, false otherwise
%
% -------
% OUTPUT:
% -------
%   x       - (1×1 double) root of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all    - (1×(k+1) double) root estimates at all iterations
%       • a_all    - (1×(k+1) double) bracketing interval lower bounds at 
%                    all iterations
%       • b_all    - (1×(k+1) double) bracketing interval upper bounds at
%                    at all iterations
%       • nintiter - (1×1 double) number of iterations to find a bracketing
%                    interval
%       • niter    - (1×1 double) number of solver iterations
%       • nfeval   - (1×1 double) number of function evaluations
%
%==========================================================================
function [x,output] = root_bisection(f,x0,opts)
    
    % sets tolerance (defaults to 10⁻¹⁰)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'TOL')
        TOL = 1e-10;
    else
        TOL = opts.TOL;
    end
    
    % sets maximum number of iterations (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'maxiter')
        maxiter = 200;
    else
        maxiter = opts.k_max;
    end
    
    % determines if the initial bracketing interval should be updated
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'rebracket')
        rebracket = false;
    else
        rebracket = opts.rebracket;
    end
    
    % starts counting function evaluations and number of iterations for
    % finding an initial bracketing interval
    nfeval = 0;
    nintiter = 0;
    
    % obtains bracketing interval
    if length(x0) == 1
        [a,b,nfeval_bracket_1,nintiter_1] = bracket_sign_change(f,x0);
        nfeval = nfeval + nfeval_bracket_1;
        nintiter = nintiter + nintiter_1;
        rebracket = false;
    else
        a = x0(1);
        b = x0(2);
    end
    
    % updates bracketing interval if requested
    if rebracket
        [a,b,nfeval_bracket_2,nintiter_2] = bracket_sign_change(f,[a,b]);
        nfeval = nfeval + nfeval_bracket_2;
        nintiter = nintiter + nintiter_2;
    end
    
    % initial guess for root
    c = (a+b)/2;
    
    % evaluates function at initial guess
    fc = f(c);
    nfeval = nfeval + 1;
    
    % returns root estimate at first iteration if it is a root of f(x)
    if fc == 0
        x = c;
        output.x_all = x;
        output.a_all = a;
        output.b_all = b;
        output.nintiter = nintiter;
        output.niter = 0;
        output.nfeval = nfeval;
        return
    end
    
    % evaluates function at lower bound of initial bracketing interval
    fa = f(a);
    
    % preallocates arrays to store all intermediate solutions and
    % bracketing intervals
    x_all = zeros(1,maxiter+1);
    a_all = zeros(1,maxiter+1);
    b_all = zeros(1,maxiter+1);
    
    % stores initial guess and bracketing interval
    x_all(1) = c;
    a_all(1) = a;
    b_all(1) = b;
    
    % prints header for solver progress
    if print
        print_solver_header;
    end
    
    % iteration
    for k = 1:maxiter
        
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
        
        % stores updated root estimate and bracketing interval
        x_all(k+1) = c;
        a_all(k+1) = a;
        b_all(k+1) = b;
        
        % terminates if converged
        if ((b-a) < TOL)
            if print, fprintf(''); end
            break;
        end
        
        % evaluates function at updated root estimate
        fc = f(c);
        nfeval = nfeval+1;
        
        % prints solver progress
        if print
            print_solver_progress(k,nfeval,c,f)
        end
        
    end
    
    % converged root
    x = c;
    
    % additional outputs
    output.x_all = x_all(1:(k+1));
    output.a_all = a_all(1:(k+1));
    output.b_all = b_all(1:(k+1));
    output.nintiter = nintiter;
    output.niter = k;
    output.nfeval = nfeval;
    
end