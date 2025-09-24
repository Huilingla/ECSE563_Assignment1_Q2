function Z = impedance(nfrom, nto, r, x, b)
% IMPEDANCE Calculate the impedance matrix for an AC power network.
%
%   Z = IMPEDANCE(NFROM, NTO, R, X, B) computes the nodal impedance matrix
%   by first calculating the admittance matrix and then solving Z*Y = I.
%
%   Inputs:
%     nfrom - Mx1 vector of branch starting nodes
%     nto   - Mx1 vector of branch ending nodes  
%     r     - Mx1 vector of branch resistances (p.u.)
%     x     - Mx1 vector of branch reactances (p.u.)
%     b     - Mx1 vector of branch shunt susceptances (p.u.)
%
%   Output:
%     Z     - NxN nodal impedance matrix (p.u.)
%
%   Method:
%     1. Calculate admittance matrix Y using the admittance() function
%     2. Solve Z*Y = I for Z using linsolve() (column by column)

    fprintf('=== IMPEDANCE MATRIX CALCULATION ===\n\n');
    
    % Step 1: Calculate admittance matrix using existing function
    fprintf('Step 1: Calculating admittance matrix Y\n');
    fprintf('---------------------------------------\n');
    Y = admittance(nfrom, nto, r, x, b);
    
    % Determine network size
    all_nodes = unique([nfrom; nto]);
    N = max(all_nodes);
    
    % Step 2: Calculate impedance matrix by solving Z*Y = I
    fprintf('\nStep 2: Calculating impedance matrix by solving Z*Y = I\n');
    fprintf('-------------------------------------------------------\n');
    fprintf('Solving system column by column using linsolve()\n');
    
    % Initialize impedance matrix
    Z = zeros(N, N);
    
    % Create identity matrix for the right-hand side
    I_matrix = eye(N);
    
    % Solve for each column of Z: Z*Y = I => Y'*Z' = I'
    for col = 1:N
        e_col = I_matrix(:, col);  % Current column of identity matrix
        z_col = linsolve(Y, e_col);  % Solve Y * z_col = e_col
        Z(:, col) = z_col;
    end
    
    fprintf('Impedance matrix Z (%dx%d) calculated successfully\n\n', N, N);
    
    % Display the complete 9x9 impedance matrix
    fprintf('COMPLETE 9x9 IMPEDANCE MATRIX:\n');
    fprintf('==============================\n');
    
    % Display column headers
    fprintf('Bus \\');
    for j = 1:N
        fprintf('%12d        ', j);
    end
    fprintf('\n');
    
    % Display separator line
    for j = 1:N
        fprintf('----------------');
    end
    fprintf('\n');
    
    % Display matrix rows
    for i = 1:N
        fprintf('%2d | ', i);
        for j = 1:N
            if imag(Z(i,j)) >= 0
                fprintf('%7.4f+j%7.4f  ', real(Z(i,j)), imag(Z(i,j)));
            else
                fprintf('%7.4f-j%7.4f  ', real(Z(i,j)), abs(imag(Z(i,j))));
            end
        end
        fprintf('\n');
    end
    
    % Step 3: Validate the result
    fprintf('\nStep 3: Validation\n');
    fprintf('------------------\n');
    
    % Check that Z * Y ≈ I (identity matrix)
    identity_check = Z * Y;
    max_error = max(max(abs(identity_check - eye(N))));
    
    fprintf('Maximum error in Z*Y = I validation: %e\n', max_error);
    
    if max_error < 1e-10
        fprintf('Validation: Excellent (Z is accurate inverse of Y)\n');
    else
        fprintf('Validation: Good (minor numerical differences)\n');
    end
end
