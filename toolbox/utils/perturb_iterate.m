%==========================================================================
%
% perturb_iterate  Perturb an iterate.
%
%   xp = perturb_iterate(x)
%   xp = perturb_iterate(x,h)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-19
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
%   x       - (n×1 double) iterate
%   h       - (OPTIONAL) (1×1 double) relative step size (defaults to 100ε)
%
% -------
% OUTPUT:
% -------
%   xp      - (n×1 double) perturbed iterate
%
%==========================================================================
function xp = perturb_iterate(x,h)
    
    % defaults relative step size to 100ε if not input
    if (nargin < 2) || isempty(h)
        h = 100*eps;
    end
    
    % perturbs iterate
    xp = x+h*(1+abs(x));
    
end