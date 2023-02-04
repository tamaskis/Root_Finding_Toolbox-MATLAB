%==========================================================================
%
% print_solver_header  Prints the header for the progress printout of an
% iterative solver.
%
%   print_solver_header
%
% Copyright © 2022 Tamas Kis
% Last Update: 2023-02-04
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/Root_Finding_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Root_Finding_Methods.pdf
%
%==========================================================================
function print_solver_header
    fprintf('\n%9s     %10s     %11s     %11s\n','Iteration',...
        'Func-count','x','f(x)');
end