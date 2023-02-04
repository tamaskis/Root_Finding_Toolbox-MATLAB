%==========================================================================
%
% convergence_plot  Generates a convergence plot to investigate the 
% convergence of various finite difference methods in solving the one-
% dimensional heat equation.
%
%   convergence_plot(physical,conditions,solver,convergence)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-07-04
%
%--------------------------------------------------------------------------
%
% MATLAB Central File Exchange: 
% GitHub: https://github.com/tamaskis/heat_equation_1D-MATLAB
%
% See EXAMPLES.mlx for examples and "DOCUMENTATION.pdf" for additional 
% documentation. Both of these files are included with the download.
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   physical    - (struct) physical parameters -- fields:
%                   D   - (1×1) diffusivity
%                   a   - (1×1) left domain boundary
%                   b   - (1×1) right domain boundary
%   conditions  - (struct) initial/boundary conditions -- fields:
%                   f0      - (function_handle) f0(x) = initial condition
%                   g       - (function_handle) g(t) = left boundary cond.
%                   g_type  - (char) 'Dirichlet' or 'Neumann'
%                   h       - (function_handle) h(t) = right boundary cond.
%                   h_type  - (char) 'Dirichlet' or 'Neumann'
%   solver      - (struct) solver parameters -- fields:
%                   method  - (char) 'FTCS', 'FTCS simultaneous', 
%                             'FTCS looped', 'BTCS', 'Crank-Nicolson', or
%                             'CTCS'
%                   dt      - (1×1) time step
%                   tf      - (1×1) simulation end time
%                   dx      - (1×1) grid spacing
%                   N       - (1×1) number of subintervals
%   convergence - (struct) convergence plot parameters -- fields:
%                   Nr      - (1×1) number of refinement iterations
%                   N0      - (1×1) initial number of subintervals
%                   dt0     - (1×1) initial time step
%
% -----
% NOTE:
% -----
%   --> The "convergence" structure should NOT have both the "N0" and "dt0"
%       fields. If "N0" is supplied, "convergence_plot" examines
%       convergence with respect to the grid spacing. If "dt0" is supplied,
%       "convergence_plot" examines convergence with respect to the time
%       step.
%   --> Nt = length of time vector
%   --> N = number of subintervals (for mesh)
%
%==========================================================================
function convergence_plot(physical,conditions,solver,convergence)
    
    % extracts "a" and "b"
    a = physical.a;
    b = physical.b;
    
    % booleans that indicate whether we are refining the mesh or time step
    refining_mesh = false;
    refining_time_step = false;
    
    % extracts parameters for convergence analysis
    Nr = convergence.Nr;
    if isfield(convergence,'N0')
        N0 = convergence.N0;
        refining_mesh = true;
    elseif isfield(convergence,'dt0')
        dt0 = convergence.dt0;
        refining_time_step = true;
    end
    
    % preallocates arrays
    p = zeros(Nr-1,1);
    e = zeros(Nr,1);
    f = zeros(Nr+1,1);
    N = zeros(Nr+1,1);
    dx = zeros(Nr+1,1);
    dt = zeros(Nr+1,1);
    
    % initialize parameters at first refinement iteration
    if refining_mesh
        N(1) = N0;
        dx(1) = (b-a)/N(1);
    elseif refining_time_step
        N(1) = solver.N;
        dt(1) = dt0;
    end
    
    % calculates solutions for each refinement iteration
    for i = 1:(Nr+1)
        
        % performs refinement
        if i > 1
            if refining_mesh
                N(i) = N(i-1)*2;
                dx(i) = (b-a)/N(i);
            elseif refining_time_step
                N(i) = N(i-1);
                dt(i) = dt(i-1)/2;
            end
        end
        
        % stores relevant parameter in "solver" structure
        if refining_mesh
            solver.N = N(i);
        elseif refining_time_step
            solver.dt = dt(i);
        end

        % calculates solution
        fi = heat_equation_1D(physical,conditions,solver);

        % stores solution at central node at final time
        f(i) = fi(end,N(i)/2+1);
        
    end
    
    % calculates approximate error (i.e. error between iterations)
    for i = 1:Nr
        e(i) = f(i)-f(i+1);
    end
    
    % calculates convergence rate
    for i = 1:(Nr-1)
        p(i) = log2(e(i)/e(i+1));
    end
    
    % determines which convergence rate it is closer to
    if abs(mean(p)-2) < abs(mean(p)-1)
        p_true = 2;
    else
        p_true = 1;
    end
    
    % parameters for convergence plot
    if refining_mesh
        x = dx;
        xlabel_str = 'Mesh Spacing, $\Delta x$';
    elseif refining_time_step
        x = dt;
        xlabel_str = 'Time Step, $\Delta t$';
    end
    
    % scaling term for plotting
    if p_true == 1
        de = x(end)/abs(e(end));
    else
        de = x(end)^2/abs(e(end));
    end
    
    % convergence plot
    figure('Position',[475,350,600,500]);
    loglog(x,x,'k:','linewidth',1.5);
    xlim([x(Nr),x(1)]);
    hold on;
    loglog(x,x.^2,'k--','linewidth',1.5);
    loglog(x(1:Nr),abs(e)*de,'linewidth',1.5,'color','k');
    grid on;
    xlabel(xlabel_str,'interpreter','latex','fontsize',18);
    ylabel("Approximate Error Magnitude, $|\varepsilon|$",'interpreter',...
        'latex','fontsize',18);
    legend("$1^{\textrm{st}}\textrm{ Order Convergence}$",...
        "$2^{\textrm{nd}}\textrm{ Order Convergence}$",...
        "$\textrm{Calculated Convergence}$",'interpreter','latex',...
        'fontsize',14,'location','south');
    hold off;

end