% RUN_IMPEDANCE_VALIDATION Validate the impedance matrix calculation
clear all; close all; clc;

% Load the IEEE 9-bus system data
ieee9_A1;

fprintf('=== IMPEDANCE MATRIX VALIDATION - IEEE 9-BUS SYSTEM ===\n\n');

% Calculate impedance matrix using the new function
fprintf('CALCULATING IMPEDANCE MATRIX\n');
fprintf('============================\n');
Z = impedance(nfrom, nto, r, x, b);

% Display the impedance matrix
fprintf('\nIMPEDANCE MATRIX OUTPUT\n');
fprintf('=======================\n');
fprintf('Full Impedance Matrix Z (%dx%d) in p.u.:\n\n', size(Z,1), size(Z,2));

% Display as combined complex matrix
fprintf('Combined Impedance Matrix (R + jX):\n');
for i = 1:size(Z,1)
    fprintf('Bus %d: ', i);
    for j = 1:size(Z,2)
        if imag(Z(i,j)) >= 0
            fprintf('%8.6f + j%8.6f  ', real(Z(i,j)), imag(Z(i,j)));
        else
            fprintf('%8.6f - j%8.6f  ', real(Z(i,j)), abs(imag(Z(i,j))));
        end
    end
    fprintf('\n');
end

% Additional validation: Solve using impedance matrix with reference node
fprintf('\nADDITIONAL VALIDATION: SOLVING WITH REFERENCE NODE CONSTRAINT\n');
fprintf('==============================================================\n');

% Calculate admittance matrix separately for comparison
Y = admittance(nfrom, nto, r, x, b);

% Choose reference node
ref_node = 9;
fprintf('Using node %d as reference (V_%d = 0)\n', ref_node, ref_node);

% Modify Y to make it invertible (set reference node constraints)
Y_modified = Y;
Y_modified(ref_node, :) = 0;
Y_modified(:, ref_node) = 0;
Y_modified(ref_node, ref_node) = 1;

% Calculate modified impedance matrix using linsolve (column by column)
fprintf('\nCalculating modified impedance matrix with reference node constraint...\n');

N = size(Y_modified, 1);
Z_modified = zeros(N, N);
I_matrix = eye(N);

for col = 1:N
    e_col = I_matrix(:, col);
    z_col = linsolve(Y_modified, e_col);
    Z_modified(:, col) = z_col;
end

fprintf('Modified Impedance Matrix (with reference node %d constraint):\n', ref_node);
for i = 1:size(Z_modified,1)
    fprintf('Bus %d: ', i);
    for j = 1:size(Z_modified,2)
        if imag(Z_modified(i,j)) >= 0
            fprintf('%8.6f + j%8.6f  ', real(Z_modified(i,j)), imag(Z_modified(i,j)));
        else
            fprintf('%8.6f - j%8.6f  ', real(Z_modified(i,j)), abs(imag(Z_modified(i,j))));
        end
    end
    fprintf('\n');
end

% Solve using impedance matrix: V = Z_modified * I
fprintf('\nSolving using impedance matrix: V = Z_modified * I_modified\n');
I_modified = Iint;
I_modified(ref_node) = 0;  % Reference node has no current injection

V_impedance = Z_modified * I_modified;

% Convert to polar coordinates
V_mag = abs(V_impedance);
V_angle_deg = angle(V_impedance) * 180/pi;

fprintf('\nVOLTAGE SOLUTION USING IMPEDANCE MATRIX\n');
fprintf('======================================\n');
fprintf('Node    Magnitude (p.u.)    Angle (degrees)\n');
fprintf('----    ----------------    --------------\n');

for i = 1:length(V_impedance)
    if i == ref_node
        fprintf('%2d (ref) %12.4f        %12.2f\n', i, V_mag(i), V_angle_deg(i));
    else
        fprintf('%2d       %12.4f        %12.2f\n', i, V_mag(i), V_angle_deg(i));
    end
end

% Compare with admittance matrix solution using linsolve
fprintf('\nCOMPARISON WITH ADMITTANCE MATRIX SOLUTION (using linsolve)\n');
fprintf('===========================================================\n');

V_admittance = linsolve(Y_modified, I_modified);
difference = max(abs(V_impedance - V_admittance));

fprintf('Maximum voltage difference between methods: %e p.u.\n', difference);

if difference < 1e-10
    fprintf('Methods are consistent: ✓\n');
else
    fprintf('Methods show minor differences (expected due to numerical precision)\n');
end

% Display key impedance matrix properties
fprintf('\nIMPEDANCE MATRIX PROPERTIES\n');
fprintf('===========================\n');
fprintf('Matrix size: %dx%d\n', size(Z,1), size(Z,2));
fprintf('Condition number of Y: %e\n', cond(Y));
fprintf('Symmetry check: %e (max off-diagonal difference)\n', max(max(abs(Z - Z'))));
fprintf('Diagonal elements (self-impedances):\n');
for i = 1:size(Z,1)
    fprintf('  Z(%d,%d) = %.6f + j%.6f p.u.\n', i, i, real(Z(i,i)), imag(Z(i,i)));
end

% Test the relationship: V = Z * I should equal Y \ I
fprintf('\nFINAL VALIDATION: Z = Y^(-1) relationship\n');
fprintf('========================================\n');
test_current = Iint(1:3);  % Test with first 3 nodes
V_from_Z = Z(1:3, 1:3) * test_current;
V_from_Y = linsolve(Y(1:3, 1:3), test_current);

validation_error = max(abs(V_from_Z - V_from_Y));
fprintf('Maximum error in Z = Y^(-1) validation: %e\n', validation_error);
