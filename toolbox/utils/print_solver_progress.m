%==========================================================================
%
% print_solver_progress  Print the progress of an iterative solver.
%
%   print_solver_progress(iter,nfeval,x)
%   print_solver_progress(iter,nfeval,x,f)
%   print_solver_progress(iter,nfeval,x,f,a,b)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2023-02-25
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
%   iter    - (1×1 double) iteration number
%   nfeval  - (1×1 double) number of function evaluations
%   x       - (1×1 double) current iterate
%   f       - (OPTIONAL) (1×1 double) function evaluation at current
%             iterate
%   a       - (OPTIONAL) (1×1 double) lower bound of current interval
%   b       - (OPTIONAL) (1×1 double) upper bound of current interval
%
%==========================================================================
function print_solver_progress(iter,nfeval,x,f,a,b)
    
    % determines if function evaluations and intervals are provided
    f_input = (nargin >= 4) && ~isempty(f);
    a_input = (nargin >= 5) && ~isempty(a);
    b_input = (nargin == 6) && ~isempty(b);
    
    % format for iteration number
    iter_format = '%9i     ';
    
    % format for number of function evaluations
    nfeval_format = '%10i     ';
    
    % format for lower bound of current interval
    if a_input
        if a < 0
            a_format = '%1.4e     ';
        else
            a_format = ' %1.4e     ';
        end
    end
    
    % format for upper bound of current interval
    if b_input
        if b < 0
            b_format = '%1.4e     ';
        else
            b_format = ' %1.4e     ';
        end
    end
    
    % format for current iterate
    if x < 0
        x_format = '%1.4e     ';
    else
        x_format = ' %1.4e     ';
    end
    if ~f_input
        x_format = [x_format,'\n'];
    end
    
    % format for function evaluation at current iterate
    if f < 0
        f_format = '%1.4e     \n';
    else
        f_format = ' %1.4e\n';
    end
    
    % prints solver progress
    if f_input
        if a_input && b_input
            progress_format = [iter_format,nfeval_format,a_format,...
                b_format,x_format,f_format];
            fprintf(progress_format,iter,nfeval,a,b,x,f);
        else
            progress_format = [iter_format,nfeval_format,x_format,...
                f_format];
            fprintf(progress_format,iter,nfeval,x,f);
        end
    else
        if a_input && b_input
            progress_format = [iter_format,nfeval_format,a_format,...
                b_format,x_format];
            fprintf(progress_format,iter,nfeval,a,b,x);
        else
            progress_format = [iter_format,nfeval_format,x_format];
            fprintf(progress_format,iter,nfeval,x);
        end
    end
    
end