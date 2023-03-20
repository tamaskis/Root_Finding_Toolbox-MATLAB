%==========================================================================
%
% rootn_broyden  Broyden's method for finding the root of a multivariate, 
% vector-valued function.
%
%   x = rootn_broyden(f,J,x0)
%   x = rootn_broyden(f,J,x0,opts)
%   [x,output] = rootn_broyden(__)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-01-08
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/Root_Finding_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Root_Finding_Methods.pdf
%
% DEPENDENCIES:
%   • Numerical Differentiation Toolbox (https://tamaskis.github.io/Numerical_Differentiation_Toolbox-MATLAB/)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) multivariate, vector-valued function,
%             f(x) (f : ℝⁿ → ℝⁿ)
%   x0      - (n×1 double) initial guess for root
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • J          - (1×1 function_handle) Jacobian of f(x)
%                      (J : ℝⁿ → ℝⁿˣⁿ)
%       • xatol      - (1×1 double) absolute step tolerance (defaults to 
%                      10⁻¹⁰)
%       • vtol       - (1×1 double) value tolerance (defaults to 0)
%       • max_iter   - (1×1 double) maximimum number solver of iterations
%                      allowed (defaults to 200)
%       • max_feval  - (1×1 double) maximimum number of function
%                      evaluations allowed (defaults to 200n TODO)
%       • print      - (1×1 logical) true if solver progress should be
%                      printed, false otherwise (defaults to false)
%
% -------
% OUTPUT:
% -------
%   x       - (n×1 double) root of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all      - (n×(n_iter+1) double) root estimates at all 
%                      iterations
%       • f_all      - (n×(n_iter+1) double) function evaluations at all
%                      iterations
%       • n_iter     - (1×1 double) number of solver iterations
%       • n_feval    - (1×1 double) number of function evaluations
%       • n_jeval    - (1×1 double) number of Jacobian evaluations
%
% -----
% NOTE:
% -----
%   --> n_feval does NOT count any function evaluations performed
%       internally by a user-specified Jacobian (via opts.J).
%
%==========================================================================
function [x,output] = rootn_broyden(f,x0,opts)
    
    % TODO: shouldn't need 200n, just 200
    % dimension of x
    n = length(x0);
    
    % sets absolute step tolerance (defaults to 10⁻¹⁰)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'xatol')
        xatol = 1e-10;
    else
        xatol = opts.xatol;
    end
    
    % sets value tolerance (defaults to 0ₙ)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'vtol')
        vtol = zeros(n,1);
    else
        vtol = opts.vtol;
    end
    
    % sets maximum number of iterations allowed (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'max_iter')
        max_iter = 200;
    else
        max_iter = opts.max_iter;
    end
    
    % sets maximum number of function evaluations allowed (defaults to 
    % 200n) TODO
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'max_feval')
        max_feval = 200*n;
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
    
    % initializes previous and next root estimates
    x_prev = x0;
    x_next = 0;
    
    % evaluates function at initial guess
    f_prev = f(x_prev);
    n_feval = n_feval+1;
    
    % returns initial guess if it is a root of f(x)
    if f_prev <= vtol
        x = x_prev;
        output.x_all = x_prev;
        output.f_all = f_prev;
        output.n_iter = 0;
        output.n_feval = n_feval;
        return
    end
    
    % evaluates Jacobian of function at initial guess (defaults to using 
    % central difference approximation)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'J')
        J0 = cjacobian(f,x0);
        n_feval = n_feval+2*n;
    else
        J0 = opts.J(x0);
    end
    
    % inverse of the initial Jacobian
    A = J0\eye(n);
    
    % Newton step
    s = A*f_prev;
    
    % root estimate and function evaluation at first iteration
    x_curr = x0+s;
    f_curr = f(x_curr);
    n_feval = n_feval+1;
    
    % preallocates arrays to store root estimates and function evaluations
    x_all = zeros(n,max_iter+1);
    f_all = zeros(n,max_iter+1);
    
    % stores initial guess and first iterate and their corresponding
    % function evaluations
    x_all(:,1) = x_prev;
    x_all(:,2) = x_curr;
    f_all(:,1) = f_prev;
    f_all(:,2) = f_curr;
    
    % prints header for solver progress
    if print
        print_solver_header;
    end
    
    % iteration
    for k = 2:max_iter
        
        % updates the inverse Jacobian via the Sherman-Morrison formula
        y = f_curr-f_prev;
        z = -A*y;
        p = -s.'*z;
        uT = s.'*A;
        A = A+((s+z)*uT)/p;
        
        % Newton step
        s = -A*f_curr;
        
        % updates root estimate
        x_next = x_curr+s;
        
        % updates function evaluation for next iteration
        f_next = f(x_next);
        n_feval = n_feval+1;
        
        % stores kth root estimate and corresponding function evaluation
        x_all(:,k+1) = x_next;
        f_all(:,k+1) = f_next;
        
%         % prints solver progress TODO
%         if print
%             print_solver_progress(k,n_feval,x_next,f_next)
%         end
        
        % solver termination
        if all(abs(f_next) < vtol) || (norm(s) <= xatol) ||...
                (n_feval >= max_feval)
            break;
        end
        
        % stores updated values for next iteration
        x_curr = x_next;
        f_prev = f_curr;
        f_curr = f_next;
        
    end
    
    % converged root
    x = x_next;
    
    % number of iterations
    n_iter = k;
    
    % additional outputs
    output.x_all = x_all(:,1:(n_iter+1));
    output.f_all = f_all(:,1:(n_iter+1));
    output.n_iter = n_iter;
    output.n_feval = n_feval;
    
end