%==========================================================================
%
% rootn_iteration  Function iteration method for finding the root of a 
% multivariate, vector-valued function.
%
%   x = rootn_iteration(f,x0)
%   x = rootn_iteration(f,x0,opts)
%   [x,k] = rootn_iteration(__)
%   [x,k,x_all] = rootn_iteration(__)
%
% See also TODO.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-12-18
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
%   x0      - (1×1 double) initial guess for root
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
%   x       - (1×1 double) root of f(x)
%   k       - (1×1 double) number of iterations
%   x_all   - (1×(k+1) double) root estimates at all iterations
%
%==========================================================================
function [x,k,x_all] = rootn_iteration(f,x0,opts)
    
    % sets tolerance (defaults to 10⁻¹⁰)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'TOL')
        opts.TOL = 1e-10;
    end
    
    % sets maximum number of iterations (defaults to 200)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'k_max')
        opts.k_max = 200;
    end
    
    % determines if all intermediate estimates should be returned
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'return_all')
        opts.return_all = false;
    end
    
    % defines the auxiliary function
    g = @(x) x-f(x);
    
    % solves for the root of f(x) using fixed-point iteration on g(x)
    if opts.return_all
        [x,k,x_all] = fixed_point_n(g,x0,opts);
    else
        [x,k] = fixed_point_n(g,x0,opts);
    end
    
end