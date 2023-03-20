%==========================================================================
%
% root_newton  Newton's method for finding the root of a differentiable,
% univariate, scalar-valued function.
%
%   x = root_newton(f,df,x0)
%   x = root_newton(f,df,x0,opts)
%   [x,output] = root_newton(__)
%
% See also root_bisection, root_brent_dekker, root_iteration, root_itp,
% root_secant.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-18
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
%   df      - (1×1 function_handle) derivative of f(x) (f' : ℝ → ℝ)
%   x0      - (1×1 double) initial guess for root
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • xatol      - (1×1 double) absolute step tolerance (defaults to 
%                      10⁻¹⁰)
%       • vtol       - (1×1 double) value tolerance (defaults to 0)
%       • max_iter   - (1×1 double) maximimum number solver of iterations
%                      allowed (defaults to 200)
%       • max_feval  - (1×1 double) maximimum number of function
%                      evaluations allowed (defaults to 200)
%       • max_deval  - (1×1 double) maximimum number of derivative
%                      evaluations allowed (defaults to 200)
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
%       • f_all      - (1×(n_iter+1) double) function evaluations at all
%                      iterations
%       • df_all     - (1×(n_iter+1) double) derivative evaluations at all
%                      iterations
%       • n_iter     - (1×1 double) number of solver iterations
%       • n_feval    - (1×1 double) number of function evaluations
%       • n_deval    - (1×1 double) number of derivative evaluations
%
%==========================================================================
function [x,output] = root_newton(f,df,x0,opts)
    
    % sets absolute step tolerance (defaults to 10⁻¹⁰)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'xatol')
        xatol = 1e-10;
    else
        xatol = opts.xatol;
    end
    
    % sets value tolerance (defaults to 0)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'vtol')
        vtol = 0;
    else
        vtol = opts.vtol;
    end
    
    % sets maximum number of iterations allowed (defaults to 200)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'max_iter')
        max_iter = 200;
    else
        max_iter = opts.max_iter;
    end
    
    % sets maximum number of function evaluations allowed (defaults to 200)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'max_feval')
        max_feval = 200;
    else
        max_feval = opts.max_feval;
    end
    
    % sets maximum number of derivative evaluations allowed (defaults to 
    % 200)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'max_deval')
        max_deval = 200;
    else
        max_deval = opts.max_deval;
    end
    
    % determines if solver progress should be printed (defaults to no)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'print')
        print = false;
    else
        print = opts.print;
    end
    
    % starts counting function and derivative evaluations
    n_feval = 0;
    n_deval = 0;
    
    % initializes current and next root estimates
    x_curr = x0;
    x_next = 0;
    
    % evaluates function and derivative at initial guess
    f_curr = f(x_curr);
    df_curr = df(x_curr);
    n_feval = n_feval+1;
    n_deval = n_deval+1;
    
    % returns initial guess if it is a root of f(x)
    if abs(f_curr) <= vtol
        x = x_curr;
        output.x_all = x_curr;
        output.f_all = f_curr;
        output.df_all = df_curr;
        output.n_iter = 0;
        output.n_feval = n_feval;
        output.n_deval = n_deval;
        return
    end
    
    % preallocates arrays to store root estimates, function evaluations,
    % and derivative evaluations
    x_all = zeros(1,max_iter+1);
    f_all = zeros(1,max_iter+1);
    df_all = zeros(1,max_iter+1);
    
    % stores initial guess and its corresponding function and derivative
    % evaluations
    x_all(1) = x_curr;
    f_all(1) = f_curr;
    df_all(1) = df_curr;
    
    % prints header for solver progress
    if print
        print_solver_header;
    end
    
    % iteration
    for k = 1:max_iter
        
        % perturbs current root estimate if derivative is 0 and 
        % re-evaluates function and its derivative at perturbed root 
        % estimate
        if df_curr == 0
            x_curr = perturb_iterate(x_curr);
            f_curr = f(x_curr);
            df_curr = df(x_curr);
            n_feval = n_feval+1;
            n_deval = n_deval+1;
        end
        
        % updates root estimate
        x_next = x_curr-f_curr/df_curr;
        
        % evaluates function and its derivative at updated root estimate
        f_next = f(x_next);
        df_next = df(x_next);
        n_feval = n_feval+1;
        n_deval = n_deval+1;
        
        % stores kth root estimate and its corresponding function and
        % derivative evaluations
        x_all(k+1) = x_next;
        f_all(k+1) = f_next;
        df_all(k+1) = df_next;
        
        % prints solver progress
        if print
            print_solver_progress(k,n_feval,x_next,f_next)
        end
        
        % solver termination on convergence criteria
        if (abs(f_next) < vtol) || (abs(x_next-x_curr) <= xatol)
            break;
        end
        
        % solver termination on timeout criteria
        if (n_feval >= max_feval) || (n_deval >= max_deval)
            break;
        end
        
        % stores updated values for next iteration
        x_curr = x_next;
        f_curr = f_next;
        df_curr = df_next;
        
    end
    
    % converged root
    x = x_next;
    
    % number of iterations
    n_iter = k;
    
    % additional outputs
    output.x_all = x_all(1:(n_iter+1));
    output.f_all = f_all(1:(n_iter+1));
    output.df_all = df_all(1:(n_iter+1));
    output.n_iter = n_iter;
    output.n_feval = n_feval;
    output.n_deval = n_deval;
    
end