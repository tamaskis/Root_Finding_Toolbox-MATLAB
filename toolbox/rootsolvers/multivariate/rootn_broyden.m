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
% DEPENDENCIES:
%   • Numerical Differentiation Toolbox (https://www.mathworks.com/matlabcentral/fileexchange/97267-numerical-differentiation-toolbox)
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
%       • TOL        - (1×1 double) tolerance (defaults to 10⁻¹⁰)
%       • k_max      - (1×1 double) maximimum number of iterations, kₘₐₓ
%                      (defaults to 200)
%       • J          - (1×1 function_handle) Jacobian of f(x) 
%                      (J : ℝⁿ → ℝⁿˣⁿ) (defaults to using the Jacobian 
%                      approximation provided by a differentiator object)
%       • d          - (1×1 Differentiator) differentiator object (defaults
%                      to a default differentiator object)
%
% -------
% OUTPUT:
% -------
%   x       - (n×1 double) root of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all   - (n×(k+1) double) root estimates at all iterations
%       • k       - (1×1 double) number of solver iterations
%       • f_count - (1×1 double) number of function evaluations
%
%==========================================================================
function [x,output] = rootn_broyden(f,x0,opts)
    
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
    
    % determines if a Jacobian is input, and extracts it if it is
    jacobian_input = (nargin == 3) && ~isempty(opts) && isfield(opts,'J');
    if jacobian_input
        J = opts.J;
    end
    
    % determines if a differentiator object is input, and extracts it if it
    % is
    diff_input = (nargin == 3) && ~isempty(opts) && isfield(opts,'d');
    if diff_input
        d = opts.d;
    end
    
    % creates a default differentiator object if needed
    if ~jacobian_input && ~diff_input
        d = Differentiator();
    end
    
    % dimension of x
    n = length(x0);
    
    % evaluates function at initial guess
    f_prev = f(x0);
    
    % returns initial guess if it is a root of f(x)
    if f_prev == zeros(n,1)
        x = x0;
        output.x_all = x;
        output.k = 0;
        output.f_count = 1;
        return
    end
    
    % evaluates Jacobian of function at initial guess
    if jacobian_input
        J0 = J(x0);
    else
        J0 = d.jacobian(f,x0);
    end
    
    % inverse of the initial Jacobian
    %A = inv(J0); TODO
    A = J0\eye(n);
    
    % Newton step
    s = A*f_prev;
    
    % root estimate at first iteration
    x_curr = x0+s;
    
    % preallocates array to store all intermediate solutions
    x_all = zeros(n,k_max+1);
    
    % stores initial guess and root estimate at 1st iteration
    x_all(:,1) = x0;
    x_all(:,2) = x_curr;
    
    % iteration
    for k = 2:k_max
        
        % evaluates function at current iteration
        f_curr = f(x_curr);
        
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
        
        % stores updated root estimate
        x_all(:,k+1) = x_next;
        
        % terminates solver if converged
        if (norm(s) < TOL)
            break;
        end
        
        % stores updated root estimate and current function evaluation for 
        % next iteration
        x_curr = x_next;
        f_prev = f_curr;
        
    end
    
    % converged root
    x = x_next;
    
    % additional outputs
    output.x_all = x_all(:,1:(k+1));
    output.k = k;
    output.f_count = k+1;
    
end