%==========================================================================
%
% fixedpointn_iteration  Fixed-point iteration for finding the fixed point
% of a multivariate, vector-valued function.
%
%   c = fixedpointn_iteration(f,x0)
%   c = fixedpointn_iteration(f,x0,opts)
%   [c,output] = fixedpointn_iteration(__)
%
% See also fixed_point.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-01-07
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
%   x0      - (n×1 double) initial guess for fixed point
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
%   c       - (n×1 double) fixed point of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • c_all   - (n×(k+1) double) fixed point estimates at all 
%                   iterations
%       • k       - (1×1 double) number of solver iterations
%       • f_count - (1×1 double) number of function evaluations
%
%==========================================================================
function [c,output] = fixedpointn_iteration(f,x0,opts)
    
    % sets tolerance (defaults to 10⁻¹⁰)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'TOL')
        TOL = 1e-10;
    else
        TOL = opts.TOL;
    end
    
    % sets maximum number of iterations (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'k_max')
        k_max = 200;
    else
        k_max = opts.k_max;
    end
    
    % returns initial guess if it is a fixed point of f(x)
    if norm(f(x0)-x0) < TOL
        c = x0;
        output.c_all = c;
        output.k = 0;
        output.f_count = 1;
        return
    end
    
    % dimension of x
    n = length(x0);
    
    % initializes current and next fixed point estimates
    x_curr = x0;
    x_next = zeros(n,1);
    
    % preallocates array to store all intermediate solutions and stores
    % initial guess
    c_all = zeros(n,k_max+1);
    c_all(:,1) = x0;
    
    % iteration
    for k = 1:k_max
        
        % updates fixed point estimate
        x_next = f(x_curr);
        
        % stores updated fixed point estimate
        c_all(:,k+1) = x_next;
        
        % terminates if converged
        if (norm(x_next-x_curr) < TOL)
            break;
        end
        
        % stores updated fixed point estimate for next iteration
        x_curr = x_next;
        
    end
    
    % prints blank line after last line of solver progress printouts
    if print, fprintf(''); end
    
    % converged fixed point
    c = x_next;
    
    % additional outputs
    output.c_all = c_all(:,1:(k+1));
    output.k = k;
    output.f_count = k+1;
    
end