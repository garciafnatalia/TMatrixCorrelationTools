clear; 
%%%%%%%%%%%%%%%%%%%%% SYSTEM B - Reference system %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Well-separated modes %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of system B 
wn2_0=diag([200 1000 2500]);
B0 = [1 -1 0.5; 1 1 -1; -1 0.5 1]'; 
M=inv(B0*B0');                      % Mass matrix
K=inv(B0*inv(wn2_0)*B0');           % Stiffness matrix
[B,wn2_B]=eig(K,M);                 % Modal matrix B and omega^2
clear wn2_0 B0

%%%%%%%%%%%%%% SYSTEM A - Perturbed system (STIFFNESS CHANGE) %%%%%%%%%%%%%
DK = [1.0 -0.5 0.0;-0.5 2.0 0.5;0.0 0.5 1.5 ]*15;   % Stiffness change
DM=[0 0 0; 0 0 0; 0 0 0];                        % Mass change
[auxA,auxwn2_A]=eig(K+DK,M+DM);
[d,ind] = sort(diag(auxwn2_A));
wn2_A =auxwn2_A(ind,ind);                        % Omega^2 of system A
% Mass normalization of mode shapes:
A = auxA(:,ind);
for s=1:2
    A(:,s)=A(:,s)/sqrt(A(:,s)'*(M+DM)*A(:,s));
end
clear auxwn2_A auxA d ind s 

%%%%%%%%%%%%%%%%%%%%%%%%%% Correlation analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%
T=inv(B)*A;

T_Mass_Matrix = T_Mass(B,A) 
%To visualize, use: T_Mass_Matrix = T_Mass(B,A,'plot') 

F=sqrt(diag(wn2_B))/(2*pi)
[T_Stiffness_Matrix1,T_Stiffness_Matrix2] = T_Stiffness(B,A,F)
%[T_Stiffness_Matrix1,T_Stiffness_Matrix2] = T_Stiffness(B,A,F,'plot')
[MAC_Matrix,ROTMAC_Matrix] = Rotmac(B, A)
%[MAC_Matrix,ROTMAC_Matrix] = Rotmac(B, A, 'plot')

%INTERPRETATION:
%  - T_Mass_Matrix values equal to 90° indicate perfect mass correlation. 
%  - T_Stiffness_Matrix detects stiffness discrepancies. Note that 
%    T_Stiffness_Matrix1 is more sensitive to stiffness changes but also to 
%    identification errors in the mode shapes.
%  - MAC values show slight discrepancies due to rotations. ROTMAC displays 
%    diagonal values equal to 1, which confirms the absence of shear (as 
%    expected when there are only stiffness changes).
%  - Similar results are expected in the case of closely-spaced modes, 
%    since in the case of stiffness discrepancies only rotational effects 
%    appear, and shear effects are absent.

