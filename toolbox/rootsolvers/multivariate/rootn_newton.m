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
%       • TOL        - (1×1 double) tolerance (defaults to 10⁻¹⁰)
%       • k_max      - (1×1 double) maximimum number of iterations, kₘₐₓ
%                      (defaults to 200)
%
% -------
% OUTPUT:
% -------
%   x       - (n×1 double) root of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all   - (n×(k+1) double) root estimates at all iterations
%       • k       - (1×1 double) number of solver iterations
%       • f_count - (1×1 double) number of function evaluations
%       • J_count - (1×1 double) number of Jacobian evaluations
%
%==========================================================================
function [x,output] = rootn_newton(f,J,x0,opts)
    
    % sets tolerance (defaults to 10⁻¹⁰)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'TOL')
        TOL = 1e-10;
    else
        TOL = opts.TOL;
    end
    
    % sets maximum number of iterations (defaults to 200)
    if (nargin < 4) || isempty(opts) || ~isfield(opts,'k_max')
        k_max = 200;
    else
        k_max = opts.k_max;
    end
    
    % turns singular matrix warning into error
    s = warning('error','MATLAB:singularMatrix');
    
    % dimension of x
    n = length(x0);
    
    % returns initial guess if it is a root of f(x)
    if f(x0) == zeros(n,1)
        x = x0;
        output.x_all = x;
        output.k = 0;
        output.f_count = 1;
        output.J_count = 0;
        return
    end
    
    % sets solution estimate at the first iteration of Newton's method
    % as the initial guess
    x_curr = x0;
    
    % zero vector
    zero_vector = zeros(n,1);
    
    % initializes x_new so its scope isn't limited to the while loop
    x_next = zero_vector;
    
    % preallocates array to store all intermediate solutions and stores
    % initial guess
    x_all = zeros(n,k_max+1);
    x_all(:,1) = x0;
    
    % counter for number of times current root estimate is perturbed
    n_perturb = 0;
    
    % iteration
    for k = 1:k_max
        
        % solves for y assuming Jacobian is nonsingular
        try
            f_curr = f(x_curr);
            y = J(x_curr)\(-f_curr);
            
        % catch any error when solving linear system
        catch ME
            
            % perturbs the initial guess if exception was due to Jacobian
            % singularity
            if strcmpi(ME.identifier,'MATLAB:singularMatrix')
                
                % stores unperturbed current root estimate in case we need
                % to return it
                x_curr_unperturbed = x_curr;
                
                % perturbs current root estimate
                if x_curr ~= zero_vector
                    x_curr = x_curr*(1+100*TOL*abs(x_curr));
                else
                    x_curr = 100*TOL;
                end
                
                % tries to resolve system at perturbed root estimate
                try
                    n_perturb = n_perturb+1;
                    y = J(x_curr)\(-f(x_curr));
                    
                % if solution of perturbed system fails, check whether the
                % function evaluation at unperturbed system is near 0 
                % (occasionally the Jacobian will be singular at the root)
                catch ME
                    if f_curr < TOL
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
        if (norm(y) < TOL)
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