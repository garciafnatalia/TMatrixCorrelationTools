function T_Mass_Matrix = T_Mass(varargin)
% FUNCTION DESCRIPTION:
%   Computes the T-Mass indicator to assess correlation in terms of mass  
%   between two models.
%
% USAGE:
%   T_Mass_Matrix = T_Mass(M_B, M_A)
%   T_Mass_Matrix = T_Mass(M_B, M_A,'plot')
%   T_Mass_Matrix = T_Mass(T)
%   T_Mass_Matrix = T_Mass(T,'plot')
%
% INPUTS:
%   M_B, M_A  - Modal matrices of systems B and A (modes in columns).
%   T         - Directly provided transformation matrix T [T = pinv(B)*A].
%   'plot'    - (Optional) Include this keyword to generate 3 figures.
%                 
%   *Important notes
%       - It is strongly recommended to review the theoretical background 
%    	  before using these indicators.
%       - M_B must be mass-normalized.
%       - When comparing a numerical model with an experimental one, M_B 
%         should correspond to the numerical modal matrix.
%       - T can be directly provided by the user, for example by assembling  
%         it from submatrices to reduce truncation effects.
%
% OUTPUTS:
%   T_Mass_Matrix - T-Mass correlation indicator: Angles between the column 
%                   vectors of matrix T. Off-diagonal terms equal to 90º
%                   indicate perfect correlation in terms of mass.
%   FIGURES:
%       Figure(1) - Heatmap of the T-Mass matrix.
%       Figure(2) - T-Mass values [°] sorted in descending order.
%       Figure(3) - T-Mass values [°] grouped by column. Each point is
%                   labeled with its matrix index (row, column).
%
% REFERENCE:
% When using the outcome of this function in scientific publications,
% please cite:
%   [1] García Fernández, N., Fernández Fernandez, P., Brincker, R., & 
%   Aenlle López, M. (2024). Mass and Stiffness Correlation Using a 
%   Transformation Matrix. Infrastructures, 9(6), 96. 
%   https://doi.org/10.3390/infrastructures9060096 
%
%
% Implementation by NATALIA GARCÍA FERNÁNDEZ (garciafnatalia@uniovi.es)
%-------------------------------------------------------------------------

plot_flag = any(strcmpi(varargin, 'plot'));
args = varargin(~strcmpi(varargin, 'plot'));  

switch numel(args)
    case 1
        T = varargin{1};
    case 2
        B = varargin{1};
        A = varargin{2};

        [rB, cB] = size(B);
        [rA, cA] = size(A);

        if rB ~= rA
            error('Both modal matrices must have the same number of rows (DOFs).');
        end

        T = pinv(B) * A;
    otherwise
        error('Invalid input. Usage: T_Mass(M_B,M_A) or T_Mass(T) [with optional ''plot''].');
end

[~, cT] = size(T);
T_Mass_Matrix = zeros(cT, cT);

for p = 1:cT
    for q = 1:cT
        T_Mass_Matrix(p,q)=(180/pi)*subspace(T(:,p),T(:,q));
    end
end
T_Mass_Matrix(logical(eye(size(T_Mass_Matrix)))) = NaN;

if plot_flag

    %FIGURE 1
    figure()    
    vmin = min(T_Mass_Matrix(~isnan(T_Mass_Matrix)));
    vmax = max(T_Mass_Matrix(~isnan(T_Mass_Matrix)));
    if abs(vmax - vmin) < 1e-4
        if abs(vmax) ==90
        vmax = 90;
        vmin = 0;
        else
        vmin = vmin - 0.01;
        vmax = vmax + 0.01;
        end
    end    
    h = heatmap(T_Mass_Matrix);
    h.Title = 'T-Mass matrix [º]';
    h.Colormap = parula %jet
    h.ColorLimits = [vmin, vmax];             
    h.MissingDataColor = [1 1 1];
    h.MissingDataLabel = 'NaN'; 
    h.CellLabelFormat = '%.3f';

    %FIGURE 2
    mask = ~eye(size(T_Mass_Matrix));
    T_Mass_values = T_Mass_Matrix(mask);
    [y, idx] = sort(T_Mass_values, 'descend');
    x = 1:length(y);
    figure;
    plot(x, y, 'o:','Color', 'k','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    xlabel('Matrix elements (sorted)');
    ylabel('T-Mass values [°]');
    yr = ylim;
    yrange = abs(yr(2) - yr(1));
    if yrange < 0.001
        ytickformat('%.3f');
    end

    %FIGURE 3
    [n, ~] = size(T_Mass_Matrix);
    colors = lines(n);  
    x_vals = [];    
    y_vals = [];    
    labels = {};     
    col_idx = [];    
    idx = 1;
    for col = 1:n
        for row = 1:n
            if row ~= col
                x_vals(end+1) = idx;
                y_vals(end+1) = T_Mass_Matrix(row, col);
                labels{end+1} = sprintf('(%d,%d)', row, col);
                col_idx(end+1) = col;
                idx = idx + 1;
            end
        end
    end

    figure; hold on;
    for col = 1:n
        mask = col_idx == col;
        plot(x_vals(mask), y_vals(mask), 'o:', ...
             'MarkerEdgeColor', colors(col, :), ...
             'MarkerFaceColor', colors(col, :), ...
             'Color', colors(col, :), ...
             'LineWidth', 1.2);
    end
    xticks(x_vals);
    xticklabels(labels);
    xtickangle(45); 
    xlabel('Matrix index');
    ylabel('T-Mass values [°]');
    grid on;
    box on;
    yr = ylim;
    yrange = abs(yr(2) - yr(1));
    if yrange < 0.001
        ytickformat('%.3f');
    end
    
end
    
end
