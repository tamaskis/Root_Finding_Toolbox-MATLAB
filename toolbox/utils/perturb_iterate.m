%==========================================================================
%
% perturb_iterate  Perturb an iterate.
%
%   xp = perturb_iterate(x)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-03
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
%
% -------
% OUTPUT:
% -------
%   xp      - (n×1 double) perturbed iterate
%
%==========================================================================
function xp = perturb_iterate(x)
    
    % iterate dimension
    n = length(x);
    
    % perturbs iterate
    if x ~= zeros(n,1)
        xp = x*(1+max(100*eps,eps*norm(x)));
    else
        xp = 100*eps*ones(n,1);
    end
    
end