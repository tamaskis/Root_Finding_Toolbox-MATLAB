%==========================================================================
%
% root_iteration  Function iteration method for finding the root of a 
% univariate, scalar-valued function.
%
%   x = root_iteration(f,x0)
%   x = root_iteration(f,x0,opts)
%   [x,output] = root_iteration(__)
%
% See also root_bisection, root_brent_dekker, root_itp, root_newton, 
% root_secant.
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
%   f       - (1×1 function_handle) univariate, scalar-valued function, 
%             f(x) (f : ℝ → ℝ)
%   x0      - (1×1 double) initial guess for root
%   opts    - (OPTIONAL) (1×1 struct) solver options
%       • TOL        - (1×1 double) tolerance (defaults to 10⁻¹⁰)
%       • k_max      - (1×1 double) maximimum number of iterations, kₘₐₓ
%                      (defaults to 200)
%
% -------
% OUTPUT:
% -------
%   x       - (1×1 double) root of f(x)
%   output  - (1×1 struct) algorithm outputs
%       • x_all   - (1×(k+1) double) root estimates at all iterations
%       • k       - (1×1 double) number of solver iterations
%       • f_count - (1×1 double) number of function evaluations
%
%==========================================================================
function [x,output] = root_iteration(f,x0,opts)
    
    % defaults opts to empty array if not input
    if (nargin < 3)
        opts = [];
    end
    
    % defines the auxiliary function
    g = @(x) x-f(x);
    
    % solves for the root of f(x) using fixed-point iteration on g(x)
    [x,output] = fixed_point(g,x0,opts);
    
    % renames "c_all" field of output to "x_all"
    output.x_all = output.c_all;
    output = rmfield(output,'c_all');
    
end