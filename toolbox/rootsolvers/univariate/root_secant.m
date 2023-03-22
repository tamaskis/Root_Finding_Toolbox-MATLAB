%==========================================================================
%
% root_secant  Secant method for finding the root of a univariate, 
% scalar-valued function.
%
%   x = root_secant(f,x0)
%   x = root_secant(f,x0,opts)
%   [x,output] = root_secant(__)
%
% See also root_bisection, root_brent_dekker, root_iteration, root_itp,
% root_newton.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-19
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
%   x0      - (1×1 OR 1×2 double) initial guess (x₀ ∈ ℝ) or initial guess +
%             1st iterate (x0 = [x₀,x₁] ∈ ℝ¹ˣ²)
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • xatol      - (1×1 double) absolute step tolerance (defaults to 
%                      10⁻¹⁰)
%       • vtol       - (1×1 double) value tolerance (defaults to 0)
%       • max_iter   - (1×1 double) maximimum number solver of iterations
%                      allowed (defaults to 200)
%       • max_feval  - (1×1 double) maximimum number of function
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
%       • n_iter     - (1×1 double) number of solver iterations
%       • n_feval    - (1×1 double) number of function evaluations
%
%==========================================================================
function [x,output] = root_secant(f,x0,opts)
    
    % sets absolute step tolerance (defaults to 10⁻¹⁰)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'xatol')
        xatol = 1e-10;
    else
        xatol = opts.xatol;
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
    
    % determines if solver progress should be printed (defaults to no)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'print')
        print = false;
    else
        print = opts.print;
    end
    
    % starts counting function evaluations
    n_feval = 0;
    
    % initializes previous and current root estimates
    if length(x0) == 2
        x_prev = x0(1);
        x_curr = x0(2);
    else
        x_prev = x0;
        x_curr = x_prev+sqrt(eps)*(1+abs(x0));
    end
    
    % evaluates function at initial guess
    f_prev = f(x_prev);
    n_feval = n_feval+1;
    
    % returns initial guess if it is a root of f(x)
    if abs(f_prev) <= vtol
        x = x_prev;
        output.x_all = x_prev;
        output.f_all = f_prev;
        output.n_iter = 0;
        output.n_feval = n_feval;
        return
    end
    
    % evaluates function at first iterate
    f_curr = f(x_curr);
    n_feval = n_feval+1;
    
    % returns first iterate if it is a root of f(x)
    if abs(f_prev) <= vtol
        x = x_curr;
        output.x_all = [x_prev,x_curr];
        output.f_all = [f_prev,f_curr];
        output.n_iter = 1;
        output.n_feval = n_feval;
        return
    end
    
    % preallocates arrays to store root estimates and function evaluations
    x_all = zeros(1,max_iter+1);
    f_all = zeros(1,max_iter+1);
    
    % stores initial guess and first iterate and their corresponding
    % function evaluations
    x_all(1) = x_prev;
    x_all(2) = x_curr;
    f_all(1) = f_prev;
    f_all(2) = f_curr;
    
    % prints header for solver progress
    if print
        print_solver_header;
    end
    
    % iteration
    for k = 2:max_iter
        
        % perturbs current root estimate if function evaluation is same as
        % at previous root estimate and re-evaluates function at perturbed
        % root estimate
        if f_curr == f_prev
            x_curr = pertub_iterate(x_curr);
            f_curr = f(x_curr);
            n_feval = n_feval+1;
        end
        
        % updates root estimate
        x_next = (x_prev*f_curr-x_curr*f_prev)/(f_curr-f_prev);
        
        % evaluates function at updated root estimate
        f_next = f(x_next);
        n_feval = n_feval+1;
        
        % stores kth root estimate and the corresponding function
        % evaluation
        x_all(k+1) = x_next;
        f_all(k+1) = f_next;
        
        % prints solver progress
        if print
            print_solver_progress(k,n_feval,x_next,f_next)
        end
        
        % solver termination
        if (abs(f_next) < vtol) || (abs(x_next-x_curr) <= xatol) ||...
                (n_feval >= max_feval)
            break;
        end
        
        % stores updated values for next iteration
        x_prev = x_curr;
        x_curr = x_next;
        f_prev = f_curr;
        f_curr = f_next;
        
    end
    
    % prints blank line after last line of solver progress printouts
    if print, fprintf(''); end
    
    % converged root
    x = x_next;
    
    % number of iterations
    n_iter = k;
    
    % additional outputs
    output.x_all = x_all(1:(n_iter+1));
    output.f_all = f_all(1:(n_iter+1));
    output.n_iter = n_iter;
    output.n_feval = n_feval;
    
end