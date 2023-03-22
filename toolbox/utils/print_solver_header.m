%==========================================================================
%
% print_solver_header  Prints the header for the progress printout of an
% iterative solver.
%
%   print_solver_header
%   print_solver_header(bracketing,func)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2023-03-21
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
%   bracketing  - (OPTIONAL) (1×1 logical) true if printing solver progress
%                 for a bracketing method, false otherwise (defaults to 
%                 false)
%   func        - (OPTIONAL) (1×1 logical) true if including function
%                 evaluations in the solver progress, false otherwise
%                 (defaults to true)
%
%==========================================================================
function print_solver_header(bracketing,func)
    
    % defaults bracketing to false if not input
    if (nargin < 1) || isempty(bracketing)
        bracketing = false;
    end
    
    % defaults func to true if not input
    if (nargin < 1) || isempty(func)
        func = true;
    end
    
    % prints header for solver progress printouts
    if func
        if bracketing
            format = ['\n%9s     %10s     %11s     %11s     %11s     ',...
                '%11s\n'];
            fprintf(format,'Iteration','Func-count','a','b','x','f(x)');
        else
            format = '\n%9s     %10s     %11s     %11s\n';
            fprintf(format,'Iteration','Func-count','x','f(x)');
        end
    else
        if bracketing
            format = '\n%9s     %10s     %11s     %11s     %11s\n';
            fprintf(format,'Iteration','Func-count','a','b','x');
        else
            format = '\n%9s     %10s     %11s\n';
            fprintf(format,'Iteration','Func-count','x');
        end
    end
    
end