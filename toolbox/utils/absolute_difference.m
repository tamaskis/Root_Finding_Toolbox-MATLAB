%==========================================================================
%
% absolute_difference  Absolute difference between two values.
%
%   Delta_ab = absolute_difference(a,b)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-01-17
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
%   a           - (n×1 double) value
%   b           - (n×1 double) value
%
% -------
% OUTPUT:
% -------
%   Delta_ab    - (n×1 double) absolute difference between a and b
%
%==========================================================================
function Delta_ab = absolute_difference(a,b)
    Delta_ab = norm(b-a);
end