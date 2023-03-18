%==========================================================================
%
% convergence_plot  Convergence plot.
%
%   convergence_plot(x_all)
%   fig = convergence_plot(x_all)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-05
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
%   x_all   - (n×(n_iter+1) double) all iterates
%
% TODO: stuff with fig
%==========================================================================
function fig = convergence_plot(x_all)
    
    % -------------
    % Calculations.
    % -------------
    
    % dimension of iterates
    n = size(x_all,1);
    
    % number of iterations
    n_iter = size(x_all,2)-1;
    
    % preallocates array to store normed distance between iterates
    y = zeros(1,n_iter+1);
    
    % initialize distance at first iteration to NaN
    y(1) = NaN;
    
    % normed distance between each pair of iterates
    for k = 2:(n_iter+1)
        y(k) = norm(x_all(k)-x_all(k-1));
    end
    
    % iteration numbers
    k = 1:(n_iter+1);
    
    % ---------
    % Plotting.
    % ---------
    
    % initializes new figure
    figure('position',[540,300,700,500]);
    
    % creates plot
    plot(k,y,'-o','LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor',...
        'white')
    
    % turns grid on
    grid on;
    
    % sets vertical axis to logarithmic scale
    set(gca,'YScale','log')
    
    % axis labels
    xlabel('$k$','Interpreter','latex','FontSize',18);
    if n > 1
        ylabel('$\left\|\mathbf{x}_{k}-\mathbf{x}_{k-1}\right\|$',...
            'Interpreter','latex','FontSize',18);
    else
        ylabel('$\left|x_{k}-x_{k-1}\right|$','Interpreter','latex',...
            'FontSize',18);
    end
    
end