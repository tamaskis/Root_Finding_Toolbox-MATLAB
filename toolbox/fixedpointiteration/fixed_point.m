%==========================================================================
%
% fixed_point  Fixed-point iteration for finding the fixed point of a 
% univariate, scalar-valued function.
%
%   c = fixed_point(f,x0)
%   c = fixed_point(f,x0,opts)
%   [c,k] = fixed_point(__)
%   [c,k,c_all] = fixed_point(__)
%
% See also fixed_point_n.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-12-11
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
%       • TOL        - (1×1 double) tolerance (defaults to 10⁻¹⁰)
%       • k_max      - (1×1 double) maximimum number of iterations, kₘₐₓ
%                      (defaults to 200)
%       • return_all - (1×1 logical) returns estimates at all iterations if
%                      set to "true" (defaults to false)
%
% -------
% OUTPUT:
% -------
%   c       - (1×1 double) fixed point of f(x)
%   k       - (1×1 double) number of iterations
%   c_all   - (1×(k+1) double) fixed point estimates at all iterations
%
%==========================================================================
function [c,k,c_all] = fixed_point(f,x0,opts)
    
    % sets maximum number of iterations (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'k_max')
        k_max = 200;
    else
        k_max = opts.k_max;
    end
    
    % determines if all intermediate estimates should be returned
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'return_all')
        return_all = false;
    else
        return_all = opts.return_all;
    end
    
    % sets tolerance (defaults to 10⁻¹⁰)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'TOL')
        TOL = 1e-10;
    else
        TOL = opts.TOL;
    end
    
    % returns initial guess if it is a fixed point of f(x)
    if abs(f(x0)-x0) < TOL
        c = x0;
        return
    end
    
    % fixed point estimate at first iteration
    x_curr = x0;
    
    % preallocates array to store all intermediate solutions
    if return_all
        c_all = zeros(1,k_max+1);
    end
    
    % iteration
    for k = 1:k_max
        
        % stores results
        if return_all
            c_all(k) = x_curr;
        end
        
        % updates fixed point estimate
        x_next = f(x_curr);
        
        % terminates if converged
        if (abs(x_next-x_curr) < TOL)
            break;
        end
        
        % stores updated fixed point estimate for next iteration
        x_curr = x_next;
        
    end
    
    % converged fixed point
    c = x_next;
    
    % stores converged result and trims array
    if return_all
        c_all(k+1) = c;
        c_all = c_all(1:(k+1));
    end
    
end