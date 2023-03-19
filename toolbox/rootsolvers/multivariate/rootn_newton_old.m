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
%       • vtol       - (1×1 double) value tolerance (defaults to 0)
%       • max_iter   - (1×1 double) maximimum number solver of iterations
%                      allowed (defaults to 200)
%       • max_feval  - (1×1 double) maximimum number of function
%                      evaluations allowed (defaults to 200n TODO)
%       • max_deval  - (1×1 double) maximimum number of derivative
%                      evaluations allowed (defaults to 200)
%       • print      - (1×1 logical) true if solver progress should be
%                      printed, false otherwise (defaults to false)
%
% -------
% OUTPUT:
% -------
%   x       - (n×1 double) root of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all      - (1×(n_iter+1) double) root estimates at all 
%                      iterations
%       • f_all      - (1×(n_iter+1) double) function evaluations at all
%                      iterations
%       • J_all      - (n×n×(n_iter+1) double) Jacobian evaluations at all
%                      iterations
%       • n_iter     - (1×1 double) number of solver iterations
%       • n_feval    - (1×1 double) number of function evaluations
%       • n_jeval    - (1×1 double) number of Jacobian evaluations
%
%==========================================================================
function [x,output] = rootn_newton(f,J,x0,opts)
    
    % turns singular matrix warning into error
    s = warning('error','MATLAB:singularMatrix');
    
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
    % 200n)
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
    x_next = 0;
    
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
    
    % preallocates arrays to store iterates, function evaluations, and
    % Jacobian evaluations
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
        
        % solves for y assuming Jacobian is nonsingular
        try
            y = J_curr\(-f_curr);
            
        % catch any error when solving linear system
        catch ME
            
            % perturbs the initial guess if exception was due to Jacobian
            % singularity
            if strcmpi(ME.identifier,'MATLAB:singularMatrix')
                
                % stores unperturbed current root estimate in case we need
                % to return it
                x_curr_unperturbed = x_curr;
                
                % perturbs current root estimate
                x_curr = perturb_iterate(x_curr);
                
                % evaluates function and Jacobian at perturbed root
                % estimate
                f_curr = f(x_curr);
                J_curr = J(x_curr);
                n_feval = n_feval+1;
                n_jeval = n_jeval+1;
                
                % tries to solve system at perturbed root estimate
                try
                    y = J_curr\(-f_curr);
                    
                % if solution of perturbed system fails, check whether the
                % function evaluation at unperturbed system is near 0 
                % (occasionally the Jacobian will be singular at the root)
                catch ME
                    if f_curr <= vtol
                        x_next = x_curr_unperturbed;
                        break;
                    else
                        rethrow(ME);
                    end
                    
                end
                
            else
                rethrow(ME);
            end
        end
        
        % updates root estimate
        x_next = x_curr+y;
        
        % stores updated root estimate
        x_all(:,k+1) = x_next;
        
        % terminates solver if converged
        if (norm(y) < xatol)
            break;
        end
        
        % stores updated root estimate for next iteration
        x_curr = x_next;
        
    end
    
    % converged root
    x = x_next;
    
    % additional outputs
    output.x_all = x_all(:,1:(k+1));
    output.k = k;
    output.f_count = k+n_perturb+1;
    output.J_count = k+n_perturb;
    
    % resets singular matrix error to warning
    warning(s);
    
end