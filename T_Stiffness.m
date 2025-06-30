function [T_Stiffness_Matrix1,T_Stiffness_Matrix2] = T_Stiffness(varargin)
% FUNCTION DESCRIPTION:
%   Computes the T-Stiffness indicator to assess correlation in terms of  
%   stiffness between two models.
%
% USAGE:
%   [T_Stiffness_Matrix1,T_Stiffness_Matrix2] = T_Stiffness(M_B,M_A,F)
%   [T_Stiffness_Matrix1,T_Stiffness_Matrix2] = T_Stiffness(M_B,M_A,F,'plot')
%   [T_Stiffness_Matrix1,T_Stiffness_Matrix2] = T_Stiffness(T,F)
%   [T_Stiffness_Matrix1,T_Stiffness_Matrix2] = T_Stiffness(T,F,'plot')
%
% INPUTS:
%   M_B, M_A - Modal matrices of systems B and A (modes in columns).
%   T        - Directly provided transformation matrix T [T = pinv(B)*A].
%   F        - Column vector of natural frequencies [Hz] corresponding to 
%              model B. 
%   'plot'   - (Optional) Include this keyword to generate 6 figures.
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
%   T_Stiffness_Matrix1 - T-Stiffness correlation indicator (V1): Angles
%                         between the columns of matrix T and Omega_B^2*T.  
%                         Out of the diagonal terms equal to 90° indicate  
%                         perfect correlation in terms of stiffness.
%   T_Stiffness_Matrix2 - T-Stiffness correlation indicator (V2): Angles
%                         between the columns of matrix Omega_B*T.
%                         Out of the diagonal terms equal to 90° indicate  
%                         perfect correlation in terms of stiffness.
%                       
%   FIGURES:
%       Figure(1 and 2) - Heatmap of the T-Stiffness matrix V1 and V2.
%       Figure(3 and 4) - T-Stiffness (V1 & V2) values [°] sorted in 
%                         descending order.
%       Figure(5 and 6) - T-Stiffness (V1 & V2) values [°] grouped by column. 
%                         Each point is labeled with its matrix index (row, 
%                         column).
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
    case 2
        T = args{1};
        F = args{2};

    case 3
        B = args{1};
        A = args{2};
        F = args{3};

        [rB, cB] = size(B);
        [rA, cA] = size(A);
        if rB ~= rA
            error('M_B and M_A must have the same number of rows (DOFs).');
        end

        T = pinv(B) * A;

    otherwise
        error(['Invalid input. Usage:\n' ...
               'T_Stiffness(T, F)\n' ...
               'T_Stiffness(M_B, M_A, F)\n' ...
               'Optionally add ''plot'' at the end to enable visualization.']);                      
        
end

[~, cT] = size(T);
F=F(1:cT,1);
w_B=(diag(F)*2*pi);
w2_B=(diag(F)*2*pi).^2;

T_Stiffness_Matrix1 = zeros(cT, cT);
Temp=w2_B*T;
for p=1:cT
    for q=1:cT
        T_Stiffness_Matrix1(p,q)=(180/pi)*subspace(T(:,p),Temp(:,q));
    end
end
T_Stiffness_Matrix1(logical(eye(size(T_Stiffness_Matrix1)))) = NaN;


T_Stiffness_Matrix2 = zeros(cT, cT);
Temp=w_B*T;
for p = 1:cT
    for q = 1:cT
        T_Stiffness_Matrix2(p,q)=(180/pi)*subspace(Temp(:,p),Temp(:,q));
    end
end
T_Stiffness_Matrix2(logical(eye(size(T_Stiffness_Matrix2)))) = NaN;


if plot_flag        
    %FIGURE 1 & 2   
    figure()
    vmin = min(T_Stiffness_Matrix1(~isnan(T_Stiffness_Matrix1)));
    vmax = max(T_Stiffness_Matrix1(~isnan(T_Stiffness_Matrix1)));
    if abs(vmax - vmin) < 1e-4
        if abs(vmax) ==90
        vmax = 90;
        vmin = 0;
        else
        vmin = vmin - 0.01;
        vmax = vmax + 0.01;
        end
    end  
    h = heatmap(T_Stiffness_Matrix1);
    h.Title = 'T-Stiffness matrix V1 [º]';
    h.Colormap = parula %jet
    h.ColorLimits = [vmin, vmax]; 
    h.MissingDataColor = [1 1 1];
    h.MissingDataLabel = 'NaN';    
    h.CellLabelFormat = '%.3f';

    figure()
    vmin = min(T_Stiffness_Matrix2(~isnan(T_Stiffness_Matrix2)));
    vmax = max(T_Stiffness_Matrix2(~isnan(T_Stiffness_Matrix2)));
    if abs(vmax - vmin) < 1e-4
        if abs(vmax) ==90
        vmax = 90;
        vmin = 0;
        else
        vmin = vmin - 0.01;
        vmax = vmax + 0.01;
        end
    end  
    h = heatmap(T_Stiffness_Matrix2);
    h.Title = 'T-Stiffness matrix V2 [º]';
    h.Colormap = parula %jet
    h.ColorLimits = [vmin, vmax]; 
    h.MissingDataColor = [1 1 1];
    h.MissingDataLabel = 'NaN';
    h.CellLabelFormat = '%.3f';

    %FIGURE 3 & 4
    mask = ~eye(size(T_Stiffness_Matrix1));
    T_Stiffness_values = T_Stiffness_Matrix1(mask);
    [y, idx] = sort(T_Stiffness_values, 'descend');
    x = 1:length(y);
    figure;
    plot(x, y, 'o:','Color', 'k','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    xlabel('Matrix elements (sorted)');
    ylabel('T-Stiffness V1 values [°]');
    yr = ylim;
    yrange = abs(yr(2) - yr(1));
    if yrange < 0.001
        ytickformat('%.3f');
    end

    mask = ~eye(size(T_Stiffness_Matrix2));
    T_Stiffness_values = T_Stiffness_Matrix2(mask);
    [y, idx] = sort(T_Stiffness_values, 'descend');
    x = 1:length(y);
    figure;
    plot(x, y, 'o:','Color', 'k','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    xlabel('Matrix elements (sorted)');
    ylabel('T-Stiffness V2 values [°]');
    yr = ylim;
    yrange = abs(yr(2) - yr(1));
    if yrange < 0.001
        ytickformat('%.3f');
    end

    %FIGURE 5 & 6
    [n, ~] = size(T_Stiffness_Matrix1);
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
                y_vals(end+1) = T_Stiffness_Matrix1(row, col);
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
    ylabel('T-Stiffness V1 values [°]');
    grid on;
    box on;    
    yr = ylim;
    yrange = abs(yr(2) - yr(1));
    if yrange < 0.001
        ytickformat('%.3f');
    end

    %%%
    [n, ~] = size(T_Stiffness_Matrix2);
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
                y_vals(end+1) = T_Stiffness_Matrix2(row, col);
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
    ylabel('T-Stiffness V2 values [°]');
    grid on;
    box on;   
    yr = ylim;
    yrange = abs(yr(2) - yr(1));
    if yrange < 0.001
        ytickformat('%.3f');
    end
end    
              
end
