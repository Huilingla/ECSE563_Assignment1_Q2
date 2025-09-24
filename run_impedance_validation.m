% RUN_IMPEDANCE_VALIDATION Display the complete 9x9 impedance matrix
clear all; close all; clc;

% Load the IEEE 9-bus system data
ieee9_A1;

fprintf('=== COMPLETE 9x9 IMPEDANCE MATRIX - IEEE 9-BUS SYSTEM ===\n\n');

% Calculate and display the impedance matrix
Z = impedance(nfrom, nto, r, x, b);

% Additional detailed analysis
fprintf('\nDETAILED IMPEDANCE MATRIX ANALYSIS\n');
fprintf('==================================\n');

% Display in a more compact format
fprintf('\nCompact View (R + jX format):\n');
fprintf('=============================\n');

for i = 1:9
    fprintf('Row %d: ', i);
    for j = 1:9
        if imag(Z(i,j)) >= 0
            fprintf('(%6.4f+j%6.4f) ', real(Z(i,j)), imag(Z(i,j)));
        else
            fprintf('(%6.4f-j%6.4f) ', real(Z(i,j)), abs(imag(Z(i,j))));
        end
    end
    fprintf('\n');
end

% Display matrix properties
fprintf('\nMATRIX PROPERTIES:\n');
fprintf('==================\n');
fprintf('Matrix size: %dx%d\n', size(Z,1), size(Z,2));
fprintf('Condition number: %.2e\n', cond(Z));

% Check symmetry
symmetry_error = max(max(abs(Z - Z')));
fprintf('Symmetry error: %.2e\n', symmetry_error);

% Display diagonal elements (self-impedances)
fprintf('\nDIAGONAL ELEMENTS (SELF-IMPEDANCES):\n');
fprintf('====================================\n');
for i = 1:9
    if imag(Z(i,i)) >= 0
        fprintf('Z(%d,%d) = %8.6f + j%8.6f p.u.\n', i, i, real(Z(i,i)), imag(Z(i,i)));
    else
        fprintf('Z(%d,%d) = %8.6f - j%8.6f p.u.\n', i, i, real(Z(i,i)), abs(imag(Z(i,i))));
    end
end

% Display some key off-diagonal elements (mutual impedances)
fprintf('\nKEY OFF-DIAGONAL ELEMENTS (MUTUAL IMPEDANCES):\n');
fprintf('=============================================\n');
fprintf('Z(1,4) = %8.6f + j%8.6f p.u.\n', real(Z(1,4)), imag(Z(1,4)));
fprintf('Z(4,1) = %8.6f + j%8.6f p.u.\n', real(Z(4,1)), imag(Z(4,1)));
fprintf('Z(4,5) = %8.6f + j%8.6f p.u.\n', real(Z(4,5)), imag(Z(4,5)));
fprintf('Z(8,9) = %8.6f + j%8.6f p.u.\n', real(Z(8,9)), imag(Z(8,9)));

% Verify symmetry for a few pairs
fprintf('\nSYMMETRY VERIFICATION:\n');
fprintf('======================\n');
fprintf('Z(1,4) - Z(4,1) = %.2e\n', abs(Z(1,4) - Z(4,1)));
fprintf('Z(4,5) - Z(5,4) = %.2e\n', abs(Z(4,5) - Z(5,4)));
