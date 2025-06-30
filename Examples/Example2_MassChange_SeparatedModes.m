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

%%%%%%%%%%%%%%%% SYSTEM A - Perturbed system (MASS CHANGE) %%%%%%%%%%%%%%%%
DK=[0 0 0; 0 0 0; 0 0 0];                             % Stiffness change
DM=[0.05 0.01 0.00; 0.01 0.08 -0.02; 0.00 -0.02 0.06];% Mass change
[auxA,auxwn2_A]=eig(K+DK,M+DM);
[d,ind] = sort(diag(auxwn2_A));
wn2_A =auxwn2_A(ind,ind);                             % Omega^2 of system A
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
[MAC_Matrix,ROTMAC_Matrix] = Rotmac(B, A,'plot')

%INTERPRETATION:
%  - T_Mass_Matrix detects mass discrepancies.
%  - T_Stiffness_Matrix values close to 90° indicate no stiffness 
%    discrepancies.
%  - ROTMAC shows a slight improvement over the MAC. However,in the case of
%    mass discrepancies and well-separated modes, shear affects may arise, 
%    leading to ROTMAC diagonal values ≠ 1.

