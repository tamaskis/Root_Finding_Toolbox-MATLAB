%==========================================================================
%
% rootn_newton  Newton's method for finding the root of a differentiable,
% multivariate, vector-valued function.
%
%   x = rootn_newton(f,J,x0)
%   x = rootn_newton(f,J,x0,opts)
%   [x,output] = rootn_newton(__)
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
%   f       - (1×1 function_handle) multivariate, vector-valued function,
%             f(x) (f : ℝⁿ → ℝⁿ)
%   J       - (1×1 function_handle) Jacobian of f(x) (J : ℝⁿ → ℝⁿˣⁿ)
%   x0      - (n×1 double) initial guess for root
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • xatol      - (1×1 double) absolute step tolerance (defaults to 
%                      10⁻¹⁰)
%       • vtol       - (n×1 double) value tolerance (defaults to 0ₙ)
%       • max_iter   - (1×1 double) maximimum number solver of iterations
%                      allowed (defaults to 200)
%       • max_feval  - (1×1 double) maximimum number of function
%                      evaluations allowed (defaults to 200n TODO)
%       • max_jeval  - (1×1 double) maximimum number of Jacobian
%                      evaluations allowed (defaults to 200) TODO
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
%       • J_all      - (n×n×(n_iter+1) double) Jacobian evaluations at all
%                      iterations
%       • n_iter     - (1×1 double) number of solver iterations
%       • n_feval    - (1×1 double) number of function evaluations
%       • n_jeval    - (1×1 double) number of Jacobian evaluations
%
%==========================================================================
function [x,output] = rootn_newton(f,J,x0,opts)
    
    % TODO: shouldn't need 200n, just 200
    % dimension of x
    n = length(x0);
    
    % sets absolute step tolerance (defaults to 10⁻¹⁰)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'xatol')
        xatol = 1e-10;
    else
        xatol = opts.xatol;
    end
    
    % sets value tolerance (defaults to 0ₙ)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'vtol')
        vtol = zeros(n,1);
    else
        vtol = opts.vtol;
    end
    
    % sets maximum number of iterations allowed (defaults to 200)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'max_iter')
        max_iter = 200;
    else
        max_iter = opts.max_iter;
    end
    
    % sets maximum number of function evaluations allowed (defaults to 
    % 200n) TODO
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'max_feval')
        max_feval = 200*n;
    else
        max_feval = opts.max_feval;
    end
    
    % sets maximum number of Jacobian evaluations allowed (defaults to 200)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'max_jeval')
        max_jeval = 200;
    else
        max_jeval = opts.max_jeval;
    end
    
    % determines if solver progress should be printed (defaults to no)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'print')
        print = false;
    else
        print = opts.print;
    end
    
    % starts counting function and Jacobian evaluations
    n_feval = 0;
    n_jeval = 0;
    
    % initializes current and next root estimates
    x_curr = x0;
    x_next = zeros(n,1);
    
    % evaluates function and Jacobian at initial guess
    f_curr = f(x_curr);
    J_curr = J(x_curr);
    n_feval = n_feval+1;
    n_jeval = n_jeval+1;
    
    % returns initial guess if it is a root of f(x)
    if abs(f_curr) <= vtol
        x = x_curr;
        output.x_all = x_curr;
        output.f_all = f_curr;
        output.J_all = J_curr;
        output.n_iter = 0;
        output.n_feval = n_feval;
        output.n_jeval = n_jeval;
        return
    end
    
    % preallocates arrays to store root estimates, function evaluations,
    % and Jacobian evaluations
    x_all = zeros(n,max_iter+1);
    f_all = zeros(n,max_iter+1);
    J_all = zeros(n,n,max_iter+1);
    
    % stores initial guess, function evaluation, and Jacobian evaluation
    x_all(:,1) = x_curr;
    f_all(:,1) = f_curr;
    J_all(:,:,1) = J_curr;
    
    % prints header for solver progress
    if print
        print_solver_header;
    end
    
    % iteration
    for k = 1:max_iter
        
        % perturbs current root estimate if Jacobian is singular and 
        % re-evaluates function and its Jacobian at perturbed root estimate
        if rank(J_curr) ~= n
            x_curr = perturb_iterate(x_curr);
            f_curr = f(x_curr);
            J_curr = J(x_curr);
            n_feval = n_feval+1;
            n_jeval = n_jeval+1;
        end
        
        % solves for y
        y = J_curr\(-f_curr);
        
        % updates root estimate
        x_next = x_curr+y;
        
        % evaluates function and its Jacobian at updated root estimated
        f_next = f(x_next);
        J_next = J(x_next);
        n_feval = n_feval+1;
        n_jeval = n_jeval+1;
        
        % stores kth root estimate and its corresponding function and
        % Jacobian evaluations
        x_all(:,k+1) = x_next;
        f_all(:,k+1) = f_next;
        J_all(:,:,k+1) = J_next;
        
        % prints solver progress TODO
%         if print
%             print_solver_progress(k,n_feval,x_next,f_next)
%         end
        
        % solver termination on convergence criteria
        if all(abs(f_next) < vtol) || (norm(y) <= xatol)
            break;
        end
        
        % solver termination on timeout criteria
        if (n_feval >= max_feval) || (n_jeval >= max_jeval)
            break;
        end
        
        % stores updated values for next iteration
        x_curr = x_next;
        f_curr = f_next;
        J_curr = J_next;
        
    end
    
    % prints blank line after last line of solver progress printouts
    if print, fprintf(''); end
    
    % converged root
    x = x_next;
    
    % number of iterations
    n_iter = k;
    
    % additional outputs
    output.x_all = x_all(:,1:(n_iter+1));
    output.f_all = f_all(:,1:(n_iter+1));
    output.J_all = J_all(:,:,1:(n_iter+1));
    output.n_iter = n_iter;
    output.n_feval = n_feval;
    output.n_jeval = n_jeval;
    
end