function [MAC_Matrix,ROTMAC_Matrix] = Rotmac(M_B, M_A, varargin)
% FUNCTION DESCRIPTION:
%   Computes the MAC and ROTMAC correlation matrices between two models.
%
% USAGE:
%   [MAC_Matrix,ROTMAC_Matrix] = ROTMAC(M_B, M_A)
%   [MAC_Matrix,ROTMAC_Matrix] = ROTMAC(M_B, M_A, 'plot')
%
% INPUTS:
%   M_B, M_A - Modal matrices of systems B and A (modes in columns).
%   'plot'   - (Optional) Include this keyword to generate a heatmap for
%              the MAC and ROTMAC matrices.
% 
%   *Important notes:
%       - In this code only the QR decomposition is implemented.
%
% OUTPUTS:
%   MAC_Matrix - Modal Assurance Criterion.
%   ROTMAC_Matrix - Rotated MAC.
%   FIGURES:
%       Figure(1) - Heatmap of the MAC matrix.
%       Figure(2) - 3D bar graph of the MAC matrix.
%       Figure(3) - Heatmap of the ROTMAC matrix.
%       Figure(4) - 3D bar graph of the ROTMAC matrix.
%
% REFERENCE:
% When using the outcome of this function in scientific publications,
% please cite:
%   [1] Aenlle, M., García Fernández, N., & Fernández, P.(2024). 
%   Rotation of mode shapes in structural dynamics due to mass and 
%   stiffness perturbations. Mech Syst Signal Process, 212:111269. 
%   https://doi.org/10.1016/j.ymssp.2024.111269
%
%
% Implementation by NATALIA GARCÍA FERNÁNDEZ (garciafnatalia@uniovi.es)
%-------------------------------------------------------------------------
colorstyle=jet; %parula %hot %jet

plot_flag = any(strcmpi(varargin, 'plot'));

[rB, cB] = size(M_B);
[rA, cA] = size(M_A);
if rB ~= rA
    error('M_B and M_A must have the same number of rows (DOFs).');
end
%c = min(cB, cA);
%M_B = M_B(:, 1:c);
%M_A = M_A(:, 1:c);
MAC_Matrix = macfunc(M_B, M_A); 

T=pinv(M_B)*M_A;
[R,Q] = qr(T);
M_B_R=M_B*R;
ROTMAC_Matrix = macfunc(M_B_R, M_A); 
    
if plot_flag
    plot2D(MAC_Matrix,'MAC Matrix',colorstyle)    
    plot3D(MAC_Matrix,'MAC Matrix',colorstyle)
    
    plot2D(ROTMAC_Matrix,'ROTMAC Matrix',colorstyle)   
    plot3D(ROTMAC_Matrix,'ROTMAC Matrix',colorstyle)
    
end


end

% --- Local MAC function ---
function M = macfunc(phi1, phi2)
    M = abs(phi1' * phi2).^2 ./ ...
        (sum(abs(phi1).^2, 1)' * sum(abs(phi2).^2, 1));
end

% --- Local plot2D function ---
function plot2D(x,titleplot,colorstyle)
    figure;
    h = heatmap(x);
    h.Title = titleplot;
    h.Colormap = colorstyle; 
    %h.ColorLimits = [min(x,[],"all"), max(x,[],"all")];
    h.ColorLimits = [0, 1];
    h.CellLabelFormat = '%.3f';
end

% --- Local plot3D function ---
function plot3D(x,titleplot,colorstyle)
    figure
    colormap(colorstyle);
    b=bar3(x);
    title(titleplot, 'Interpreter', 'none');
    colorbar;
    caxis([0 1]);
    cdata_sz = size(b(1).CData) ;
    z_color = repelem(x,6,4);
    z_color = mat2cell(z_color,cdata_sz(1),ones(1,size(x,2))*cdata_sz(2));
    set(b,{'CData'},z_color.');

    [r, c] = size(x);
    for i = 1:r
        for j = 1:c
            z = x(i, j);
            % Coordinates are (x, y, z) → bars at (j, i)
            text(j, i, z + 0.02, sprintf('%.3f', z), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 8, ...
                'Color', 'k');
        end
    end
end