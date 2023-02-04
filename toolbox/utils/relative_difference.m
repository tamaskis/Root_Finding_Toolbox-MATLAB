%==========================================================================
%
% relative_difference  Relative difference between two values.
%
%   delta_ab = relative_difference(a,b)
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
%   a           - (n×1 double) reference value
%   b           - (n×1 double) value
%
% -------
% OUTPUT:
% -------
%   delta_ab    - (n×1 double) relative difference between a and b,
%                 referenced to a
%
%==========================================================================
function delta_ab = relative_difference(a,b)
    delta_ab = norm(b-a)/max(norm(a),eps);
end