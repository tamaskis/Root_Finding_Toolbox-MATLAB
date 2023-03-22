%==========================================================================
%
% fixedpoint_iteration  Fixed-point iteration for finding the fixed point 
% of a univariate, scalar-valued function.
%
%   x = fixedpoint_iteration(f,x0)
%   x = fixedpoint_iteration(f,x0,opts)
%   [x,output] = fixedpoint_iteration(__)
%
% See also fixedpointn_iteration.
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
%   x0      - (1×1 double) initial guess for fixed point
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • xatol      - (1×1 double) absolute step tolerance (defaults to 
%                      10⁻¹⁰)
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
%   x       - (1×1 double) fixed point of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all      - (1×(n_iter+1) double) fixed point estimates at all 
%                      iterations
%       • f_all      - (1×(n_iter+1) double) function evaluations at all
%                      iterations
%       • n_iter     - (1×1 double) number of solver iterations
%       • n_feval    - (1×1 double) number of function evaluations
%
%==========================================================================
function [x,output] = fixedpoint_iteration(f,x0,opts)
    
    % sets absolute step tolerance (defaults to 10⁻¹⁰)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'xatol')
        xatol = 1e-10;
    else
        xatol = opts.TOL;
    end
    
    % sets maximum number of iterations allowed (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'max_iter')
        max_iter = 200;
    else
        max_iter = opts.k_max;
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
    
    % initializes current and next fixed point estimates
    x_curr = x0;
    x_next = 0;
    
    % evaluates function at initial guess
    f_curr = f(x_curr);
    n_feval = n_feval+1;
    
    % returns initial guess if it is a fixed point of f(x)
    if abs(f_curr-x_curr) < xatol
        x = x_curr;
        output.x_all = x_curr;
        output.f_all = f_curr;
        output.n_iter = 0;
        output.n_feval = n_feval;
        return
    end
    
    % preallocates arrays to store fixed point estimate and function
    % evaluations
    x_all = zeros(1,max_iter+1);
    f_all = zeros(1,max_iter+1);
    
    % stores initial guess and its corresponding function evaluation
    x_all(1) = x_curr;
    f_all(1) = f_curr;
    
    % iteration
    for k = 1:max_iter
        
        % updates fixed point estimate
        x_next = f_curr;
        
        % evaluates function at updated fixed point estimate
        f_next = f(x_next);
        n_feval = n_feval+1;
        
        % stores kth fixed point estimate and its corresponding function 
        % evaluation
        x_all(k+1) = x_next;
        f_all(k+1) = f_next;
        
        % prints solver progress
        if print
            print_solver_progress(k,n_feval,x_next,f_next)
        end
        
        % solver termination on convergence criteria
        if (abs(x_next-x_curr) <= xatol)
            break;
        end
        
        % solver termination on timeout criteria
        if (n_feval >= max_feval)
            break;
        end
        
        % stores updated values for next iteration
        x_curr = x_next;
        f_curr = f_next;
        
    end
    
    % prints blank line after last line of solver progress printouts
    if print, fprintf(''); end
    
    % converged fixed point
    x = x_next;
    
    % number of iterations
    n_iter = k;
    
    % additional outputs
    output.x_all = x_all(1:(n_iter+1));
    output.f_all = f_all(1:(n_iter+1));
    output.n_iter = n_iter;
    output.n_feval = n_feval;
    
end